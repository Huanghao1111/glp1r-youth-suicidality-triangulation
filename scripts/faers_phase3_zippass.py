#!/usr/bin/env python3
"""FAERS phase-3 zip pass: per-quarter DRUG extraction of comparator drug classes
for the notoriety-context panel (pediatric suicidality reporting inflation check).
Classes (PS/SS role, mirroring the primary GLP-1 RA exposure definition):
- ssri:  fluoxetine/sertraline/escitalopram/citalopram/fluvoxamine/paroxetine (+brands)
- sga:   aripiprazole/risperidone/quetiapine/olanzapine/ziprasidone/lurasidone/paliperidone (+brands)
- montelukast (SINGULAIR): neuropsychiatric boxed warning 2020, pediatric use
- isotretinoin (ACCUTANE etc.): long-debated suicidality association, adolescent use
Outputs compact flagged-pid CSVs to data/faers_flags3/.
Usage: python scripts/faers_phase3_zippass.py <start_idx> <end_idx>  (manifest ok rows, 0-based)
"""
import csv, os, sys, zipfile
import duckdb

ZIP_DIR = "data/faers_zips"
OUT_DIR = "data/faers_flags3"
os.makedirs(OUT_DIR, exist_ok=True)

SSRI_PAI = "FLUOXETINE|SERTRALINE|ESCITALOPRAM|CITALOPRAM|FLUVOXAMINE|PAROXETINE"
SSRI_DNAME = (SSRI_PAI + "|PROZAC|SARAFEM|RAPIFLUX|ZOLOFT|LEXAPRO|CELEXA|LUVOX|PAXIL|BRISDELLE|PEXEVA")
SGA_PAI = "ARIPIPRAZOLE|RISPERIDONE|QUETIAPINE|OLANZAPINE|ZIPRASIDONE|LURASIDONE|PALIPERIDONE"
SGA_DNAME = (SGA_PAI + "|ABILIFY|RISPERDAL|SEROQUEL|ZYPREXA|GEODONE|LATUDA|INVEGA|VRAYLAR|REXULTI")
MONT = "MONTELUKAST|SINGULAIR"
ISO = "ISOTRETINOIN|ACCUTANE|ABSORICA|CLARAVIS|AMNESTEEM|MYORISAN|ZENATANE|EPURIS"

def member(zf, prefix):
    for n in zf.namelist():
        base = os.path.basename(n).upper()
        if base.startswith(prefix) and base.endswith(".TXT"):
            return n
    return None

def main():
    rows = []
    with open("data/faers_quarters_manifest.csv", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["status"] == "ok":
                rows.append(r)
    start, end = int(sys.argv[1]), int(sys.argv[2])
    con = duckdb.connect()
    con.execute("pragma memory_limit='6GB'")

    for r in rows[start:end]:
        quarter = r["quarter"]
        zpath = os.path.join(ZIP_DIR, os.path.basename(r["url"]))
        out = os.path.join(OUT_DIR, f"{quarter}_cmp.csv")
        if os.path.exists(out):
            print(f"{quarter}: exists, skip")
            continue
        zf = zipfile.ZipFile(zpath)
        dm = member(zf, "DRUG")
        dpath = zf.extract(dm, path="data/faers_tmp")
        with open(dpath, encoding="utf-8", errors="ignore") as f:
            header = f.readline().upper()
        has_pai = "PROD_AI" in header
        pai_expr = "upper(coalesce(prod_ai,''))" if has_pai else "''"
        sql = f"""
        copy (
          select pid, max(ssri) ssri, max(sga) sga, max(mont) mont, max(iso) iso from (
            select try_cast(primaryid as bigint) pid,
              case when upper(role_cod) in ('PS','SS') and (regexp_matches({pai_expr}, '{SSRI_PAI}') or regexp_matches(upper(coalesce(drugname,'')), '{SSRI_DNAME}')) then 1 else 0 end ssri,
              case when upper(role_cod) in ('PS','SS') and (regexp_matches({pai_expr}, '{SGA_PAI}') or regexp_matches(upper(coalesce(drugname,'')), '{SGA_DNAME}')) then 1 else 0 end sga,
              case when upper(role_cod) in ('PS','SS') and (regexp_matches({pai_expr}, 'MONTELUKAST') or regexp_matches(upper(coalesce(drugname,'')), '{MONT}')) then 1 else 0 end mont,
              case when upper(role_cod) in ('PS','SS') and (regexp_matches({pai_expr}, 'ISOTRETINOIN') or regexp_matches(upper(coalesce(drugname,'')), '{ISO}')) then 1 else 0 end iso
            from read_csv('{dpath.replace(os.sep, '/')}', delim='$', header=true, all_varchar=true, quote='', escape='')
          ) t where ssri=1 or sga=1 or mont=1 or iso=1 group by pid
        ) to '{out.replace(os.sep, '/')}' (header, delim ',')"""
        con.execute(sql)
        os.remove(dpath)
        n = sum(1 for _ in open(out)) - 1
        print(f"{quarter}: flagged pids={n}", flush=True)
    con.close()

if __name__ == "__main__":
    main()
