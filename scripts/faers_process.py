# -*- coding: utf-8 -*-
# FAERS 季度包流式处理: 解压 DEMO/DRUG/REAC/DELETE -> DuckDB 打标 -> 每季度一个 parquet
# 用法: python scripts/faers_process.py <start_idx> <end_idx>  (manifest ok 行的 0 基索引)
import csv, os, re, shutil, sys, time, zipfile
import duckdb

ZIP_DIR = "data/faers_zips"
TMP_DIR = "data/faers_tmp"
OUT_DIR = "data/faers_proc"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)

SUBSTANCES = {  # 列名 -> (prod_ai 正则, drugname 正则)
    "sema": ("SEMAGLUTIDE",  "SEMAGLUTIDE|OZEMPIC|WEGOVY|RYBELSUS"),
    "lira": ("LIRAGLUTIDE",  "LIRAGLUTIDE|VICTOZA|SAXENDA|XULTOPHY"),
    "tirz": ("TIRZEPATIDE",  "TIRZEPATIDE|MOUNJARO|ZEPBOUND"),
    "exen": ("EXENATIDE",    "EXENATIDE|BYETTA|BYDUREON"),
    "dula": ("DULAGLUTIDE",  "DULAGLUTIDE|TRULICITY"),
    "lixi": ("LIXISENATIDE", "LIXISENATIDE|ADLYXIN|LYXUMIA|SOLIQUA"),
    "setm": ("SETMELANOTIDE","SETMELANOTIDE|IMCIVREE"),
}
ANY_PAI = "|".join(v[0] for v in SUBSTANCES.values())
ANY_DNAME = "|".join(v[1] for v in SUBSTANCES.values())

EVENT_GROUPS = {
    "ev_suic": ["SUICIDAL IDEATION","SUICIDE ATTEMPT","COMPLETED SUICIDE",
                "INTENTIONAL SELF-INJURY","SELF-INJURIOUS IDEATION","SUICIDAL BEHAVIOUR"],
    "ev_depr": ["DEPRESSION","DEPRESSED MOOD","MAJOR DEPRESSION","DEPRESSIVE SYMPTOM"],
    "ev_anx":  ["ANXIETY","ANXIETY DISORDER","PANIC ATTACK","PANIC DISORDER"],
    "ev_eat":  ["ANOREXIA NERVOSA","BULIMIA NERVOSA","BINGE EATING DISORDER",
                "EATING DISORDER","AVOIDANT/RESTRICTIVE FOOD INTAKE DISORDER"],
    "ev_psy":  ["PSYCHOTIC DISORDER","HALLUCINATION","DELUSION","PSYCHOTIC SYMPTOM"],
}

def strip_bom(path):
    with open(path, "rb") as f:
        head = f.read(3)
    if head == b"\xef\xbb\xbf":
        with open(path, "rb") as f:
            data = f.read()
        with open(path, "wb") as f:
            f.write(data[3:])
        return True
    return False

def build_sql(qdir, quarter, sex_col, has_pai, outc_col):
    demo = os.path.join(qdir, "demo.txt").replace("\\", "/")
    drug = os.path.join(qdir, "drug.txt").replace("\\", "/")
    reac = os.path.join(qdir, "reac.txt").replace("\\", "/")
    indi = os.path.join(qdir, "indi.txt").replace("\\", "/")
    outc = os.path.join(qdir, "outc.txt").replace("\\", "/")
    outp = os.path.join(OUT_DIR, f"{quarter}.parquet").replace("\\", "/")

    pai_expr = "upper(coalesce(prod_ai,''))" if has_pai else "''"
    drug_cases = ",\n".join(
        f"max(case when regexp_matches(pai,'{p}') or regexp_matches(dname,'{d}') then 1 else 0 end) as {k}"
        for k, (p, d) in SUBSTANCES.items())
    ev_cols = ",\n".join(
        f"max(case when ptu in ({','.join(repr(p) for p in pts)}) then 1 else 0 end) as {k}"
        for k, pts in EVENT_GROUPS.items())
    glp_cols = ", ".join(f"coalesce(m.{k},0) as {k}" for k in SUBSTANCES)
    ev_sel  = ", ".join(f"coalesce(r.{k},0) as {k}" for k in EVENT_GROUPS)

    return f"""
create or replace temp view demo_raw as
  select * from read_csv('{demo}', delim='$', header=true, all_varchar=true, quote='', escape='');
create or replace temp view demo as
  select pid as primaryid, cid as caseid, cver as caseversion, fdt as fda_dt,
         sex_u as sex, country_u as occr_country, age_years
  from (
    select try_cast(primaryid as bigint) pid,
           try_cast(caseid as bigint) cid,
           try_cast(caseversion as integer) cver,
           try_cast(fda_dt as integer) fdt,
           upper({sex_col}) sex_u, upper(occr_country) country_u,
           case upper(age_cod)
             when 'YR'  then try_cast(age as double)
             when 'MON' then try_cast(age as double)/12.0
             when 'WK'  then try_cast(age as double)*7.0/365.25
             when 'DY'  then try_cast(age as double)/365.25
             when 'DEC' then try_cast(age as double)*10.0
             else null end as age_years
    from demo_raw
  ) t where pid is not null;

create or replace temp view drug_rows as
  select try_cast(primaryid as bigint) pid, try_cast(drug_seq as integer) dseq,
         {pai_expr} pai, upper(coalesce(drugname,'')) dname
  from read_csv('{drug}', delim='$', header=true, all_varchar=true, quote='', escape='')
  where upper(role_cod) in ('PS','SS');

create or replace temp view glp_rows as
  select * from drug_rows
  where regexp_matches(pai,'{ANY_PAI}') or regexp_matches(dname,'{ANY_DNAME}');

create or replace temp view drug_match as
  select pid, {drug_cases}
  from drug_rows
  where regexp_matches(pai,'{ANY_PAI}') or regexp_matches(dname,'{ANY_DNAME}')
  group by pid;

create or replace temp view reac_match as
with r as (
  select try_cast(primaryid as bigint) pid, upper(pt) ptu
  from read_csv('{reac}', delim='$', header=true, all_varchar=true, quote='', escape='')
)
select pid, {ev_cols}
from r group by pid;

create or replace temp view indi_match as
with i as (
  select try_cast(primaryid as bigint) pid, try_cast(indi_drug_seq as integer) dseq,
         upper(indi_pt) ipt
  from read_csv('{indi}', delim='$', header=true, all_varchar=true, quote='', escape='')
)
select g.pid, string_agg(distinct i.ipt, '|') as glp_indi
from glp_rows g join i on g.pid = i.pid and g.dseq = i.dseq
group by g.pid;

create or replace temp view outc_match as
with o as (
  select try_cast(primaryid as bigint) pid, upper({outc_col}) oc
  from read_csv('{outc}', delim='$', header=true, all_varchar=true, quote='', escape='')
)
select pid,
  max(case when oc='DE' then 1 else 0 end) as outc_de,
  max(case when oc='LT' then 1 else 0 end) as outc_lt,
  max(case when oc='HO' then 1 else 0 end) as outc_ho,
  max(case when oc='DS' then 1 else 0 end) as outc_ds,
  max(case when oc='CA' then 1 else 0 end) as outc_ca,
  max(case when oc='RI' then 1 else 0 end) as outc_ri,
  max(case when oc='OT' then 1 else 0 end) as outc_ot,
  max(case when oc in ('DE','LT','HO','DS','CA','RI') then 1 else 0 end) as serious
from o group by pid;

copy (
  select d.primaryid, d.caseid, d.caseversion, d.fda_dt, '{quarter}' as quarter,
         d.age_years, d.sex, d.occr_country,
         {glp_cols}, {ev_sel},
         coalesce(im.glp_indi, '') as glp_indi,
         coalesce(om.outc_de,0) as outc_de, coalesce(om.outc_lt,0) as outc_lt,
         coalesce(om.outc_ho,0) as outc_ho, coalesce(om.outc_ds,0) as outc_ds,
         coalesce(om.outc_ca,0) as outc_ca, coalesce(om.outc_ri,0) as outc_ri,
         coalesce(om.outc_ot,0) as outc_ot, coalesce(om.serious,0) as serious,
         case when om.pid is null then 1 else 0 end as outc_missing
  from demo d
  left join drug_match m on d.primaryid = m.pid
  left join reac_match r on d.primaryid = r.pid
  left join indi_match im on d.primaryid = im.pid
  left join outc_match om on d.primaryid = om.pid
) to '{outp}' (format parquet, compression zstd);
"""

def main():
    rows = []
    with open("data/faers_quarters_manifest.csv", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["status"] == "ok":
                rows.append(r)
    start, end = int(sys.argv[1]), int(sys.argv[2])
    con = duckdb.connect()
    con.execute("pragma memory_limit='6GB'")

    del_path = os.path.join(OUT_DIR, "deleted_caseids.txt")
    for r in rows[start:end]:
        quarter = r["quarter"]
        zip_path = os.path.join(ZIP_DIR, os.path.basename(r["url"]))
        out_parquet = os.path.join(OUT_DIR, f"{quarter}.parquet")
        if os.path.exists(out_parquet):
            print(f"跳过 {quarter} (已处理)", flush=True)
            continue
        t0 = time.time()
        qdir = os.path.join(TMP_DIR, quarter)
        os.makedirs(qdir, exist_ok=True)
        zf = zipfile.ZipFile(zip_path)
        targets = {}
        for n in zf.namelist():
            base = os.path.basename(n).upper()
            if re.match(rf"DEMO\d{{2}}Q\d(_NEW)?\.TXT$", base): targets[n] = "demo.txt"
            elif re.match(rf"DRUG\d{{2}}Q\d\.TXT$", base): targets[n] = "drug.txt"
            elif re.match(rf"REAC\d{{2}}Q\d\.TXT$", base): targets[n] = "reac.txt"
            elif re.match(rf"INDI\d{{2}}Q\d\.TXT$", base): targets[n] = "indi.txt"
            elif re.match(rf"OUTC\d{{2}}Q\d\.TXT$", base): targets[n] = "outc.txt"
            elif ("DELETEDCASES" in base or base.startswith("DELETE")) and base.endswith(".TXT"): targets[n] = "delete.txt" if "delete.txt" not in targets.values() else "delete2.txt"
        for src, dst in targets.items():
            with zf.open(src) as fi, open(os.path.join(qdir, dst), "wb") as fo:
                shutil.copyfileobj(fi, fo, 1024 * 1024)
        zf.close()
        for f_ in ("demo.txt", "drug.txt", "reac.txt", "indi.txt", "outc.txt"):
            strip_bom(os.path.join(qdir, f_))

        del_files = [os.path.join(qdir, x) for x in ("delete.txt", "delete2.txt") if os.path.exists(os.path.join(qdir, x))]
        for del_file in del_files:
            with open(del_file, encoding="latin1") as f_, \
                 open(del_path, "a", encoding="utf-8") as out_:
                for line in f_:
                    tok = line.strip().split("$")[0].strip()
                    if tok.isdigit():
                        out_.write(tok + "\n")

        with open(os.path.join(qdir, "demo.txt"), encoding="latin1") as f_:
            demo_hdr = f_.readline().strip().split("$")
        with open(os.path.join(qdir, "drug.txt"), encoding="latin1") as f_:
            drug_hdr = f_.readline().strip().split("$")
        with open(os.path.join(qdir, "outc.txt"), encoding="latin1") as f_:
            outc_hdr = f_.readline().strip().split("$")
        sex_col = "sex" if "sex" in demo_hdr else "gndr_cod"
        has_pai = "prod_ai" in drug_hdr
        outc_col = "outc_cod" if "outc_cod" in outc_hdr else "outc_code"

        con.execute(build_sql(qdir, quarter, sex_col, has_pai, outc_col))
        shutil.rmtree(qdir, ignore_errors=True)
        print(f"{quarter} 完成, 耗时 {time.time()-t0:.0f}s", flush=True)
    con.close()
    print("本批结束", flush=True)

if __name__ == "__main__":
    main()
