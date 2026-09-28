#!/usr/bin/env python3
"""FAERS phase-2 zip pass: per-quarter extraction of DRUG/REAC to derive
- glp_ps: GLP-1 RA reported as PRIMARY suspect only
- adhd:   comparator ADHD drug class (PS/SS, mirrors primary exposure definition)
- gi:     gastrointestinal control event group (NAUSEA/VOMITING/DIARRHOEA/CONSTIPATION/PANCREATITIS)
Outputs compact flagged-pid CSVs to data/faers_flags/.
Usage: python scripts/faers_phase2_zippass.py <start_idx> <end_idx>  (manifest ok rows, 0-based)
"""
import csv, os, sys, zipfile
import duckdb

ZIP_DIR = "data/faers_zips"
OUT_DIR = "data/faers_flags"
os.makedirs(OUT_DIR, exist_ok=True)

ANY_PAI = "SEMAGLUTIDE|LIRAGLUTIDE|TIRZEPATIDE|EXENATIDE|DULAGLUTIDE|LIXISENATIDE|SETMELANOTIDE"
ANY_DNAME = ("SEMAGLUTIDE|OZEMPIC|WEGOVY|RYBELSUS|LIRAGLUTIDE|VICTOZA|SAXENDA|XULTOPHY|"
             "TIRZEPATIDE|MOUNJARO|ZEPBOUND|EXENATIDE|BYETTA|BYDUREON|DULAGLUTIDE|TRULICITY|"
             "LIXISENATIDE|ADLYXIN|LYXUMIA|SOLIQUA|SETMELANOTIDE|IMCIVREE")
ADHD_PAI = "METHYLPHENIDATE|DEXMETHYLPHENIDATE|LISDEXAMFETAMINE|AMPHETAMINE|ATOMOXETINE|VILOXAZINE"
ADHD_DNAME = ("METHYLPHENIDATE|RITALIN|CONCERTA|METHYLIN|METADATE|QUILLIVANT|DAYTRANA|"
              "DEXMETHYLPHENIDATE|FOCALIN|LISDEXAMFETAMINE|VYVANSE|AMPHETAMINE|ADDERALL|"
              "DEXEDRINE|MYDAYIS|AZSTARYS|EVEKEO|ZENZEDI|DYANAVEL|ADZENYS|ATOMOXETINE|"
              "STRATTERA|VILOXAZINE|QELBREE")
GI_PTS = "('NAUSEA','VOMITING','DIARRHOEA','CONSTIPATION','PANCREATITIS')"

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
        out_drug = os.path.join(OUT_DIR, f"{quarter}_drug.csv")
        out_reac = os.path.join(OUT_DIR, f"{quarter}_gi.csv")
        if os.path.exists(out_drug) and os.path.exists(out_reac):
            print(f"{quarter}: exists, skip")
            continue
        zf = zipfile.ZipFile(zpath)
        dm, rm = member(zf, "DRUG"), member(zf, "REAC")
        dpath = zf.extract(dm, path="data/faers_tmp")
        rpath = zf.extract(rm, path="data/faers_tmp")
        # detect PROD_AI presence from header
        with open(dpath, encoding="utf-8", errors="ignore") as f:
            header = f.readline().upper()
        has_pai = "PROD_AI" in header
        pai_expr = "upper(coalesce(prod_ai,''))" if has_pai else "''"
        dsql = f"""
        copy (
          select pid, max(glp_ps) glp_ps, max(adhd) adhd from (
            select try_cast(primaryid as bigint) pid,
              case when upper(role_cod)='PS' and (regexp_matches({pai_expr}, '{ANY_PAI}') or regexp_matches(upper(coalesce(drugname,'')), '{ANY_DNAME}')) then 1 else 0 end glp_ps,
              case when upper(role_cod) in ('PS','SS') and (regexp_matches({pai_expr}, '{ADHD_PAI}') or regexp_matches(upper(coalesce(drugname,'')), '{ADHD_DNAME}')) then 1 else 0 end adhd
            from read_csv('{dpath.replace(os.sep, '/')}', delim='$', header=true, all_varchar=true, quote='', escape='')
          ) t where glp_ps=1 or adhd=1 group by pid
        ) to '{out_drug.replace(os.sep, '/')}' (header, delim ',')"""
        rsql = f"""
        copy (
          select distinct try_cast(primaryid as bigint) pid
          from read_csv('{rpath.replace(os.sep, '/')}', delim='$', header=true, all_varchar=true, quote='', escape='')
          where upper(pt) in {GI_PTS}
        ) to '{out_reac.replace(os.sep, '/')}' (header, delim ',')"""
        con.execute(dsql)
        con.execute(rsql)
        os.remove(dpath)
        os.remove(rpath)
        n1 = sum(1 for _ in open(out_drug)) - 1
        n2 = sum(1 for _ in open(out_reac)) - 1
        print(f"{quarter}: flagged pids drug={n1} gi={n2}", flush=True)
    con.close()

if __name__ == "__main__":
    main()
