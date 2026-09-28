#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Extract cis windows around GLP1R/GIPR/MC4R from child-BMI GWAS (MoBa x12, EGG)
and iPSYCH suicide-attempt GWAS (Model1/2/3). Standardized output + F stats.
hg19/GRCh37 coordinates verified via mygene.info 2026-09-27.
MC4R note: hg19 = chr18:58,038,564-58,040,001 (60.37-60.39 Mb is GRCh38!)."""
import gzip, io, os, sys
import pandas as pd
import numpy as np

SRC = "E:/CM/GLP1"
OUT = "E:/CM/GLP1/mr_input"
os.makedirs(OUT, exist_ok=True)

GENES = {
    "GLP1R": (6, 39016574, 39055519),
    "GIPR":  (19, 46171502, 46186982),
    "MC4R":  (18, 58038564, 58040001),
}
EXT = 250_000          # extraction half-window
CORE = 100_000         # flag for core +-100kb window
P_THR = 5e-8

def dist_to_gene(pos, gs, ge):
    if pos < gs: return gs - pos
    if pos > ge: return pos - ge
    return 0

def standardize(df, dataset, gene, pos_col, chr_col):
    ch, gs, ge = GENES[gene]
    lo, hi = gs - EXT, ge + EXT
    d = df[(df[chr_col] == ch) & (df[pos_col] >= lo) & (df[pos_col] <= hi)].copy()
    return d

def scan_file(path, dataset, pos_col, chr_col, colmap, sep, chunksize=2_000_000, gz=None):
    """colmap: standard->actual for rsid,ea,nea,eaf,info,beta,se,p,n"""
    hits = []
    kw = dict(sep=sep, chunksize=chunksize, dtype=str, na_values=["-nan", "nan", "NA", "."],
              on_bad_lines="skip", low_memory=False)
    if gz is True:
        kw["compression"] = "gzip"
    for chunk in pd.read_csv(path, **kw):
        # numeric coercion for filtering
        chunk["_chr"] = pd.to_numeric(chunk[chr_col], errors="coerce")
        chunk["_pos"] = pd.to_numeric(chunk[pos_col], errors="coerce")
        for gene in GENES:
            ch, gs, ge = GENES[gene]
            sub = chunk[(chunk["_chr"] == ch) & (chunk["_pos"] >= gs - EXT) & (chunk["_pos"] <= ge + EXT)]
            if len(sub):
                hits.append((gene, sub.copy()))
    rows = []
    for gene, sub in hits:
        ch, gs, ge = GENES[gene]
        out = pd.DataFrame({
            "rsid": sub[colmap["rsid"]] if colmap.get("rsid") else "",
            "chr": sub["_chr"].astype(int),
            "pos": sub["_pos"].astype(int),
            "ea": sub[colmap["ea"]].str.upper() if colmap.get("ea") else "",
            "nea": sub[colmap["nea"]].str.upper() if colmap.get("nea") else "",
            "eaf": pd.to_numeric(sub[colmap["eaf"]], errors="coerce") if colmap.get("eaf") else np.nan,
            "info": pd.to_numeric(sub[colmap["info"]], errors="coerce") if colmap.get("info") else np.nan,
            "beta": pd.to_numeric(sub[colmap["beta"]], errors="coerce"),
            "se": pd.to_numeric(sub[colmap["se"]], errors="coerce"),
            "p": pd.to_numeric(sub[colmap["p"]], errors="coerce"),
            "n": pd.to_numeric(sub[colmap["n"]], errors="coerce") if colmap.get("n") else np.nan,
        })
        out["dataset"] = dataset
        out["gene"] = gene
        out["dist_to_gene"] = [dist_to_gene(p, gs, ge) for p in out["pos"]]
        out["in_100kb"] = out["dist_to_gene"] <= CORE
        out["F"] = (out["beta"] / out["se"]) ** 2
        out["pass_p5e8"] = out["p"] < P_THR
        rows.append(out)
    return rows

JOBS = []
# MoBa x12
TPS = ["birth","6weeks","3months","6months","8months","1year","1.5years","2years","3years","5years","7years","8years"]
for tp in TPS:
    JOBS.append(dict(path=f"{SRC}/moba2022_bmi_{tp}.gz", dataset=f"MoBa2022_{tp}", gz=True,
                     sep=r"\s+", chr_col="CHR", pos_col="BP",
                     colmap=dict(rsid="RSID", ea="EA", nea="NEA", eaf="EAF", info="INFO",
                                 beta="BETA", se="SE", p="P", n="N")))
# EGG 2020
JOBS.append(dict(path=f"{SRC}/EGG2020_GCST90002409_buildGRCh37.tsv.gz", dataset="EGG2020_childBMI", gz=True,
                 sep="\t", chr_col="chromosome", pos_col="base_pair_location",
                 colmap=dict(rsid="variant_id", ea="effect_allele", nea="other_allele", eaf=None, info=None,
                             beta="beta", se="standard_error", p="p_value", n="TotalSampleSize")))
# ieu-b-40 adult BMI (OpenGWAS, GRCh37, csv from vcf)
JOBS.append(dict(path=f"{SRC}/ieu-b-40.vcf.gz.csv", dataset="AdultBMI_ieub40", gz=None,
                 sep=",", chr_col="chr", pos_col="pos",
                 colmap=dict(rsid="SNP", ea="effect allele", nea="other allele", eaf="eaf", info=None,
                             beta="beta", se="se", p="pval", n="samplesize")))
# iPSYCH Model1/2/3 (uncompressed; Model1 header differs)
JOBS.append(dict(path=f"{SRC}/Model1.txt", dataset="iPSYCH_Model1", gz=None,
                 sep="\t", chr_col="chr", pos_col="pos",
                 colmap=dict(rsid="SNP", ea="effect_allele", nea="other_allele", eaf="eaf", info="INFO",
                             beta="beta", se="se", p="pval", n=None)))
for m in ["Model2", "Model3"]:
    JOBS.append(dict(path=f"{SRC}/{m}.txt", dataset=f"iPSYCH_{m}", gz=None,
                     sep="\t", chr_col="CHR", pos_col="POS",
                     colmap=dict(rsid="SNP", ea="A1", nea="A2", eaf="FRQ", info="INFO",
                                 beta="BETA", se="SE", p="P", n=None)))

all_rows, summary = [], []
for job in JOBS:
    ds = job["dataset"]
    if all(os.path.exists(f"{OUT}/{ds}_{g}_window.csv") for g in GENES):
        print(f"skip {ds} (already done)", flush=True); continue
    if not os.path.exists(job["path"]):
        print(f"MISSING {job['path']}", flush=True); continue
    print(f"scanning {ds} ...", flush=True)
    rows = scan_file(**job)
    for out in rows:
        if len(out) == 0: continue
        gene = out["gene"].iloc[0]
        fn = f"{OUT}/{ds}_{gene}_window.csv"
        out.to_csv(fn, index=False)
        all_rows.append(out)
        core = out[out["in_100kb"]]
        sig = core[core["pass_p5e8"]]
        top = core.loc[core["p"].idxmin()] if core["p"].notna().any() else None
        summary.append(dict(dataset=ds, gene=gene,
                            n_window_250kb=len(out), n_100kb=len(core), n_p5e8_100kb=len(sig),
                            min_p_100kb=core["p"].min(),
                            top_snp=(top["rsid"] if top is not None else ""),
                            top_p=(top["p"] if top is not None else np.nan),
                            top_F=(top["F"] if top is not None else np.nan)))
    print(f"  done {ds}", flush=True)

# Rebuild summary from ALL window CSVs (incremental-safe)
import glob
summary = []
for fn in sorted(glob.glob(f"{OUT}/*_window.csv")):
    out = pd.read_csv(fn)
    if len(out) == 0: continue
    ds, gene = out["dataset"].iloc[0], out["gene"].iloc[0]
    core = out[out["in_100kb"]]
    sig = core[core["pass_p5e8"]]
    top = core.loc[core["p"].idxmin()] if core["p"].notna().any() else None
    summary.append(dict(dataset=ds, gene=gene,
                        n_window_250kb=len(out), n_100kb=len(core), n_p5e8_100kb=len(sig),
                        min_p_100kb=core["p"].min(),
                        top_snp=(top["rsid"] if top is not None else ""),
                        top_p=(top["p"] if top is not None else np.nan),
                        top_F=(top["F"] if top is not None else np.nan)))
summ = pd.DataFrame(summary)
summ.to_csv(f"{OUT}/cis_window_summary.csv", index=False)
pd.set_option("display.width", 200)
print(summ.to_string(index=False))
print("\nall outputs in", OUT)
