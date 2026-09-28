# -*- coding: utf-8 -*-
# 儿科 GLP-1 病例: 适应证分组 × 事件谱 × 严重度分析
import csv, re
import duckdb
import pandas as pd

con = duckdb.connect()
con.execute("pragma memory_limit='8GB'")
con.execute("create or replace view raw as select * from read_parquet('data/faers_proc/*.parquet')")
con.execute("""
create or replace view del as
select distinct try_cast(column0 as bigint) cid
from read_csv('data/faers_proc/deleted_caseids.txt', header=false, all_varchar=true)
where try_cast(column0 as bigint) is not null""")
con.execute("""
create or replace view cases as
with ranked as (
  select r.*, row_number() over (
    partition by caseid
    order by coalesce(fda_dt,0) desc, coalesce(caseversion,0) desc, primaryid desc) rn
  from raw r where caseid not in (select cid from del))
select * from ranked where rn=1""")

# ---------- 儿科背景严重度(对照) ----------
bg = con.execute("""
select count(*) n,
       sum(serious) serious, sum(outc_de) de, sum(outc_missing) miss
from cases where age_years<18""").fetchone()
print(f"儿科全部报告背景: N={bg[0]:,} 严重 {bg[1]:,} ({bg[1]/bg[0]*100:.1f}%) 死亡 {bg[2]:,} ({bg[2]/bg[0]*100:.2f}%) 转归缺失 {bg[3]:,}")

# ---------- 儿科 GLP-1 病例逐例取出 ----------
df = con.execute("""
select primaryid, caseid, fda_dt, age_years, sex,
       sema, lira, tirz, exen, dula, lixi, setm, glp_indi,
       ev_suic, ev_depr, ev_anx, ev_eat, ev_psy,
       serious, outc_de, outc_lt, outc_ho, outc_ds, outc_ca, outc_ri, outc_ot, outc_missing
from cases where age_years<18 and (sema+lira+tirz+exen+dula+lixi+setm)>0
""").fetchdf()
print(f"儿科 GLP-1 病例: {len(df)}")

# ---------- 适应证分组 ----------
PAT_SPECIAL  = re.compile(r"INTENTIONAL SELF-INJURY|FOETAL EXPOSURE|FETAL EXPOSURE")
PAT_SYNDROME = re.compile(r"BARDET|LAURENCE|PRADER|LEPTIN|MELANOCORTIN|POMC|PCSK|HYPOTHALAMIC OBESITY|LIPODYSTROPHY|GENE MUTATION|ALBRIGHT")
PAT_WEIGHT   = re.compile(r"WEIGHT|OVERWEIGHT|OBESITY|BINGE EATING")
PAT_DIABETES = re.compile(r"DIABET|GLUCOSE|GLYC|INSULIN")
PAT_METAB    = re.compile(r"POLYCYSTIC|HEPATIC|STEATOSIS|STEATOHEPATITIS|METABOLIC")
PAT_UNKNOWN  = re.compile(r"UNKNOWN INDICATION|OFF LABEL|UNAPPROVED|^$")

GROUPS = [("special",  "特殊暴露途径(自我伤害用药/宫内暴露)"),
          ("syndrome", "遗传/综合征性肥胖"),
          ("weight",   "体重/肥胖管理"),
          ("diabetes", "糖尿病及糖代谢"),
          ("metab",    "其他代谢(PCOS/脂肪肝等)"),
          ("other",    "其他"),
          ("unknown",  "未报告/未知")]

def classify(indi):
    text = indi or ""
    if PAT_SPECIAL.search(text):  return "special"
    if PAT_SYNDROME.search(text): return "syndrome"
    if PAT_WEIGHT.search(text):   return "weight"
    if PAT_DIABETES.search(text): return "diabetes"
    if PAT_METAB.search(text):    return "metab"
    pts = [p for p in text.split("|") if p]
    if not pts or all(PAT_UNKNOWN.search(p) for p in pts):
        return "unknown"
    return "other"

df["grp"] = df["glp_indi"].map(classify)
DRUGS = ["sema", "lira", "tirz", "exen", "dula", "lixi", "setm"]
DRUG_CN = {"sema": "Semaglutide", "lira": "Liraglutide", "tirz": "Tirzepatide",
           "exen": "Exenatide", "dula": "Dulaglutide", "lixi": "Lixisenatide",
           "setm": "Setmelanotide"}

# 组 × 药 计数
tab = df.groupby("grp").size().rename("total").reset_index()
for d in DRUGS:
    tab[DRUG_CN[d]] = tab["grp"].map(df[df[d] == 1].groupby("grp").size()).fillna(0).astype(int)
tab["group_cn"] = tab["grp"].map(dict(GROUPS))
tab = tab.set_index("grp").reindex([g for g, _ in GROUPS]).dropna(subset=["total"])
tab["total"] = tab["total"].astype(int)
cols = ["group_cn", "total"] + [DRUG_CN[d] for d in DRUGS]
tab[cols].to_csv("results/faers_full_ped_glp1_indication.csv", index=False)
print("\n=== 适应证主分组 × 药物 -> results/faers_full_ped_glp1_indication.csv")
print(tab[cols].to_string(index=False))

# ---------- 严重度: 总体 + 各事件组 ----------
rows = []
def serious_block(name, sub):
    n = len(sub)
    if n == 0:
        return None
    return [name, n, int(sub["serious"].sum()), f"{sub['serious'].mean()*100:.1f}%",
            int(sub["outc_de"].sum()), int(sub["outc_lt"].sum()), int(sub["outc_ho"].sum()),
            int(sub["outc_ds"].sum()), int(sub["outc_ri"].sum()), int(sub["outc_ot"].sum()),
            int(sub["outc_missing"].sum())]
rows.append(serious_block("GLP-1 儿科全部", df))
for ev, cn in [("ev_suic", "自杀相关"), ("ev_depr", "抑郁"), ("ev_anx", "焦虑"),
               ("ev_eat", "进食障碍"), ("ev_psy", "精神病性症状")]:
    rows.append(serious_block(f"其中: {cn}", df[df[ev] == 1]))
rows = [r for r in rows if r]
with open("results/faers_full_ped_glp1_serious.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["subset", "N", "serious_n", "serious_pct", "death", "life_threat",
                "hospitalization", "disability", "req_intervention", "other_serious_OT",
                "outcome_missing"])
    w.writerows(rows)
print("\n=== 严重度 -> results/faers_full_ped_glp1_serious.csv")
for r in rows:
    print(f"{r[0]:18s} N={r[1]:<4d} 严重 {r[2]:<3d} ({r[3]:>6}) 死亡 {r[4]:<2d} 住院 {r[6]:<3d} 危及生命 {r[5]:<2d} 转归缺失 {r[10]}")

# ---------- 适应证组 × 事件谱 ----------
rows2 = []
for g, gcn in GROUPS:
    sub = df[df["grp"] == g]
    if len(sub) == 0:
        continue
    rows2.append([gcn, len(sub),
                  int(sub["ev_suic"].sum()), int(sub["ev_depr"].sum()), int(sub["ev_anx"].sum()),
                  int(sub["ev_eat"].sum()), int(sub["ev_psy"].sum()),
                  int(sub["serious"].sum()), int(sub["outc_de"].sum())])
with open("results/faers_full_ped_glp1_group_x_event.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["indication_group", "N", "suicidality", "depression", "anxiety",
                "eating_disorder", "psychosis", "serious", "death"])
    w.writerows(rows2)
print("\n=== 适应证组 × 事件谱 -> results/faers_full_ped_glp1_group_x_event.csv")
for r in rows2:
    print(f"{r[0]:28s} N={r[1]:<4d} 自杀相关 {r[2]:<2d} 抑郁 {r[3]:<2d} 焦虑 {r[4]:<2d} 进食障碍 {r[5]:<2d} 严重 {r[7]:<3d} 死亡 {r[8]}")

# ---------- 自杀相关 20 例明细(写作用) ----------
det = df[df["ev_suic"] == 1].copy()
det["drugs"] = det.apply(lambda r: "+".join(DRUG_CN[d] for d in DRUGS if r[d] == 1), axis=1)
det["yr"] = (det["fda_dt"] // 10000).astype("Int64")
det_out = det[["primaryid", "yr", "age_years", "sex", "drugs", "grp", "glp_indi",
               "serious", "outc_de", "outc_lt", "outc_ho"]].sort_values(["yr", "primaryid"])
det_out.columns = ["primaryid", "year", "age", "sex", "drugs", "indication_group",
                   "indication_pt", "serious", "death", "life_threat", "hospitalization"]
det_out.to_csv("results/faers_full_ped_suicidality_cases.csv", index=False)
print("\n=== 儿科 GLP-1 × 自杀相关 20 例明细 -> results/faers_full_ped_suicidality_cases.csv")
print(det_out.to_string(index=False))
con.close()
