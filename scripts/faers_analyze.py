# -*- coding: utf-8 -*-
# FAERS 全量分析: 去重 -> 分层 -> ROR + IC025 -> 时间分层 -> 描述统计
import csv, math, os
import duckdb
from scipy.stats import gamma as gamma_dist

DRUGS = ["sema", "lira", "tirz", "exen", "dula", "lixi", "setm"]
DRUG_CN = {"sema": "Semaglutide", "lira": "Liraglutide", "tirz": "Tirzepatide",
           "exen": "Exenatide", "dula": "Dulaglutide", "lixi": "Lixisenatide",
           "setm": "Setmelanotide"}
EVENTS = {"ev_suic": "Suicidality", "ev_depr": "Depression", "ev_anx": "Anxiety",
          "ev_eat": "Eating disorder", "ev_psy": "Psychosis"}

con = duckdb.connect()
con.execute("pragma memory_limit='8GB'")
con.execute("create or replace view raw as select * from read_parquet('data/faers_proc/*.parquet')")
con.execute("""
create or replace view del as
select distinct try_cast(column0 as bigint) cid
from read_csv('data/faers_proc/deleted_caseids.txt', header=false, all_varchar=true)
where try_cast(column0 as bigint) is not null
""")
con.execute("""
create or replace view cases as
with ranked as (
  select r.*, row_number() over (
    partition by caseid
    order by coalesce(fda_dt,0) desc, coalesce(caseversion,0) desc, primaryid desc
  ) rn
  from raw r
  where caseid not in (select cid from del)
)
select primaryid, caseid, fda_dt, quarter, age_years, sex, occr_country,
       sema, lira, tirz, exen, dula, lixi, setm,
       case when sema+lira+tirz+exen+dula+lixi+setm>0 then 1 else 0 end as glp,
       ev_suic, ev_depr, ev_anx, ev_eat, ev_psy
from ranked where rn=1
""")

# ---------- 流程计数 ----------
flow = con.execute("""
select count(*) from raw
""").fetchone()[0]
n_deleted_rows = con.execute("select count(*) from raw where caseid in (select cid from del)").fetchone()[0]
n_cases = con.execute("select count(*) from cases").fetchone()[0]
strata = con.execute("""
select count(*) filter (age_years<18) as ped,
       count(*) filter (age_years>=18) as adult,
       count(*) filter (age_years is null) as unk
from cases
""").fetchone()
print(f"原始报告行(55季度合计): {flow:,}")
print(f"剔除DELETE清单涉及行: {n_deleted_rows:,}")
print(f"按caseid去重后病例数: {n_cases:,}")
print(f"  儿科(<18): {strata[0]:,}  成人(>=18): {strata[1]:,}  年龄缺失: {strata[2]:,}")

# ---------- ROR 主分析 ----------
def agg(where):
    sel = ["count(*) n_total", "sum(glp) n_glp"]
    for e in EVENTS:
        sel.append(f"sum({e}) n_{e}")
        sel.append(f"sum(glp*{e}) a_glp_{e}")
    for d in DRUGS:
        sel.append(f"sum({d}) n_{d}")
        for e in EVENTS:
            sel.append(f"sum({d}*{e}) a_{d}_{e}")
    q = f"select {', '.join(sel)} from cases where {where}"
    return con.execute(q).fetchdf().iloc[0]

def ror_ci(a, b, c, d):
    corr = 0.5 if min(a, b, c, d) == 0 else 0.0
    a_, b_, c_, d_ = a + corr, b + corr, c + corr, d + corr
    v = (a_ * d_) / (b_ * c_)
    se = math.sqrt(1/a_ + 1/b_ + 1/c_ + 1/d_)
    return v, math.exp(math.log(v) - 1.96*se), math.exp(math.log(v) + 1.96*se), corr > 0

def ic(a, b, c):
    # WHO-UMC 信息成分, Gamma(0.5,1) 先验, 后验 Gamma(a+0.5, 1+E)
    N_exp = (a + b) * (a + c) / max(a + b + c + 1, 1)  # 占位, 实际在外面算
    return N_exp

def ic025(a, n_drug, n_event, n_total):
    E = n_drug * n_event / n_total
    lo = math.log2(gamma_dist.ppf(0.025, a + 0.5, scale=1/(1 + E)))
    hi = math.log2(gamma_dist.ppf(0.975, a + 0.5, scale=1/(1 + E)))
    pt = math.log2((a + 0.5) / (E + 0.5))
    return pt, lo, hi

os.makedirs("results", exist_ok=True)

def run_stratum(name, where):
    r = agg(where)
    n_total = int(r["n_total"])
    out_rows = []
    targets = [("glp", "GLP1-class(all)")] + [(d, DRUG_CN[d]) for d in DRUGS]
    for dkey, dlabel in targets:
        n_drug = int(r[f"n_{dkey}"])
        for ekey, elabel in EVENTS.items():
            n_event = int(r[f"n_{ekey}"])
            a = int(r[f"a_{dkey}_{ekey}"])
            b, c = n_drug - a, n_event - a
            dd = n_total - n_drug - n_event + a
            v, lo, hi, corr = ror_ci(a, b, c, dd)
            icpt, iclo, ichi = ic025(a, n_drug, n_event, n_total)
            out_rows.append([name, dlabel, elabel, a, b, c, dd,
                             round(v, 3), round(lo, 3), round(hi, 3), corr,
                             round(icpt, 3), round(iclo, 3), round(ichi, 3),
                             lo > 1, iclo > 0])
    path = f"results/faers_full_ror_{name}.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["stratum","drug","event_group","a","b","c","d","ROR","ROR_lo","ROR_hi",
                    "haldane","IC","IC025","IC975","signal_ROR","signal_IC025"])
        w.writerows(out_rows)
    print(f"\n=== {name} (N={n_total:,}) -> {path}")
    for row in out_rows:
        if row[1] == "GLP1-class(all)" or row[3] >= 1:
            sig = "*" if row[13] else " "
            print(f"{row[1]:16s} {row[2]:15s} a={row[3]:<5d} ROR={row[7]:>7.2f} ({row[8]:.2f}-{row[9]:.2f}){sig} IC025={row[11]:>6.2f}")
    return n_total

n_ped = run_stratum("pediatric", "age_years<18")
n_adult = run_stratum("adult", "age_years>=18")

# ---------- 时间分层(儿科, 2023年前后) ----------
period_rows = []
for period, where in [("pre-2023", "age_years<18 and fda_dt<20230101"),
                      ("post-2023", "age_years<18 and fda_dt>=20230101")]:
    r = agg(where)
    n_total = int(r["n_total"]); n_drug = int(r["n_glp"])
    for ekey, elabel in EVENTS.items():
        n_event = int(r[f"n_{ekey}"]); a = int(r[f"a_glp_{ekey}"])
        b, c = n_drug - a, n_event - a
        dd = n_total - n_drug - n_event + a
        v, lo, hi, corr = ror_ci(a, b, c, dd)
        icpt, iclo, ichi = ic025(a, n_drug, n_event, n_total)
        period_rows.append([period, "GLP1-class(all)", elabel, n_total, n_drug, a, b, c, dd,
                            round(v, 3), round(lo, 3), round(hi, 3), corr,
                            round(icpt, 3), round(iclo, 3), round(ichi, 3), lo > 1])
with open("results/faers_full_ror_period.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["period","drug","event_group","N_total","N_drug","a","b","c","d",
                "ROR","ROR_lo","ROR_hi","haldane","IC","IC025","IC975","signal_ROR"])
    w.writerows(period_rows)
print("\n=== 儿科时间分层 -> results/faers_full_ror_period.csv")
for row in period_rows:
    print(f"{row[0]:9s} {row[2]:15s} N_drug={row[4]:<4d} a={row[5]:<4d} ROR={row[9]:>7.2f} ({row[10]:.2f}-{row[11]:.2f})")

# ---------- 儿科 GLP-1 报告描述 ----------
desc = con.execute("""
select floor(fda_dt/10000)::int as yr,
       sum(sema) sema, sum(lira) lira, sum(tirz) tirz, sum(exen) exen,
       sum(dula) dula, sum(lixi) lixi, sum(setm) setm,
       sum(glp) glp_any
from cases where age_years<18 and glp=1 and fda_dt is not null
group by 1 order by 1
""").fetchdf()
desc.to_csv("results/faers_full_ped_glp1_byyear.csv", index=False)
print("\n=== 儿科 GLP-1 报告数(按年, 按药) -> results/faers_full_ped_glp1_byyear.csv")
print(desc.to_string(index=False))

ped_profile = con.execute("""
select sum(sema) n_sema, sum(lira) n_lira, sum(tirz) n_tirz, sum(exen) n_exen,
       sum(dula) n_dula, sum(lixi) n_lixi, sum(setm) n_setm, sum(glp) n_any,
       round(avg(age_years) filter (glp=1), 1) avg_age,
       round(median(age_years) filter (glp=1), 1) med_age,
       sum(case when sex='F' then 1 else 0 end) filter (glp=1) n_female,
       sum(ev_suic) filter (glp=1) c_suic, sum(ev_depr) filter (glp=1) c_depr,
       sum(ev_anx) filter (glp=1) c_anx, sum(ev_eat) filter (glp=1) c_eat,
       sum(ev_psy) filter (glp=1) c_psy
from cases where age_years<18
""").fetchdf()
print("\n=== 儿科 GLP-1 报告画像")
print(ped_profile.T.to_string())
with open("results/faers_full_flow.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["item", "value"])
    w.writerow(["raw_rows_55q", flow])
    w.writerow(["rows_deleted_caseids", n_deleted_rows])
    w.writerow(["cases_after_dedup", n_cases])
    w.writerow(["pediatric_lt18", strata[0]])
    w.writerow(["adult_ge18", strata[1]])
    w.writerow(["age_unknown", strata[2]])
con.close()
print("\n流程计数 -> results/faers_full_flow.csv")
