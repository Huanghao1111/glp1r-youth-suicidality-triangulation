#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Clump candidate cis instruments (PLINK, EUR panel) and harmonize with iPSYCH outcomes.
Inputs: E:/CM/GLP1/mr_input/*_window.csv from extract_cis.py
Outputs: clumped instrument lists + harmonized MR input tables in E:/CM/GLP1/mr_input/"""
import os, subprocess, sys
import pandas as pd
import numpy as np

DIR = "E:/CM/GLP1/mr_input"
PLINK = "E:/CM/GLP1/LD/plink_win64/plink.exe"
BFILE = "E:/CM/GLP1/LD/bfile/EUR"
P_THR = 5e-8

# candidate sets: (dataset, gene, core_only)  core_only=True -> +-100kb; False -> +-250kb
CANDIDATES = [
    ("MoBa2022_3months", "GLP1R", True),
    ("MoBa2022_6months", "GLP1R", True),
    ("MoBa2022_8months", "GLP1R", True),
    ("MoBa2022_1year",   "GLP1R", True),
    ("MoBa2022_1.5years","GLP1R", True),
    ("MoBa2022_7years",  "MC4R",  True),
    ("EGG2020_childBMI", "MC4R",  False),
    ("AdultBMI_ieub40",  "GIPR",  True),
    ("AdultBMI_ieub40",  "MC4R",  True),
]

def run_clump(ds, gene, core_only):
    win = pd.read_csv(f"{DIR}/{ds}_{gene}_window.csv")
    cand = win[win["in_100kb"]] if core_only else win
    cand = cand[(cand["p"] < P_THR) & cand["rsid"].astype(str).str.startswith("rs")].copy()
    tag = f"{ds}_{gene}_{'100kb' if core_only else '250kb'}"
    if len(cand) == 0:
        print(f"{tag}: no candidate SNPs"); return None
    inp = f"{DIR}/clumpin_{tag}.txt"
    cand[["rsid", "p"]].rename(columns={"rsid": "SNP", "p": "P"}).to_csv(inp, sep=" ", index=False)
    outpref = f"{DIR}/clump_{tag}"
    cmd = [PLINK, "--bfile", BFILE, "--clump", inp, "--clump-p1", str(P_THR),
           "--clump-p2", "1", "--clump-r2", "0.3", "--clump-kb", "100",
           "--clump-snp-field", "SNP", "--clump-field", "P", "--out", outpref, "--silent"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    cl = f"{outpref}.clumped"
    if not os.path.exists(cl):
        print(f"{tag}: clump failed\n{r.stdout}\n{r.stderr}"); return None
    # fixed-width-ish: split on whitespace
    tbl = pd.read_csv(cl, sep=r"\s+", engine="python")
    idx = tbl["SNP"].tolist()
    inst = cand[cand["rsid"].isin(idx)].copy()
    inst.to_csv(f"{DIR}/instruments_{tag}.csv", index=False)
    print(f"{tag}: {len(cand)} candidates -> {len(idx)} independent instruments: {idx}")
    return inst

PALINDROMIC = {("A","T"),("T","A"),("C","G"),("G","C")}

def harmonize(inst, outcome_path, out_name):
    out = pd.read_csv(outcome_path)
    o = out.rename(columns={"beta": "beta_out", "se": "se_out", "p": "p_out", "eaf": "eaf_out"})
    o = o[["rsid", "chr", "pos", "ea", "nea", "eaf_out", "beta_out", "se_out", "p_out"]]
    # 1) rsid match
    m = inst.merge(o.drop(columns=["chr", "pos"]), on="rsid", how="inner", suffixes=("_exp", "_outc"))
    # 2) fallback: position match for SNPs with non-rs ids on either side
    hit_rsid = set(m["rsid"])
    rest = inst[~inst["rsid"].isin(hit_rsid)]
    if len(rest):
        o2 = o.rename(columns={"rsid": "rsid_out"})
        m2 = rest.merge(o2, on=["chr", "pos"], how="inner", suffixes=("_exp", "_outc"))
        if len(m2):
            m2["rsid"] = m2["rsid"] + "|" + m2["rsid_out"].astype(str)
            m2 = m2.drop(columns=["rsid_out"])
            m = pd.concat([m, m2[m.columns]], ignore_index=True)
    rows = []
    for _, r in m.iterrows():
        ea_e, nea_e = str(r["ea_exp"]).upper(), str(r["nea_exp"]).upper()   # exposure alleles
        ea_o, nea_o = str(r["ea_outc"]).upper(), str(r["nea_outc"]).upper() # outcome alleles
        if (ea_o, nea_o) == (ea_e, nea_e):
            b = r["beta_out"]; action = "aligned"
        elif (ea_o, nea_o) == (nea_e, ea_e):
            b = -r["beta_out"]; action = "flipped"
        else:
            continue  # allele mismatch
        if (ea_e, nea_e) in PALINDROMIC:
            fe, fo = r.get("eaf_exp", np.nan), r.get("eaf_out", np.nan)
            if pd.isna(fe) or pd.isna(fo) or (0.4 < fe < 0.6) or (0.4 < fo < 0.6):
                continue  # ambiguous palindromic
            if (fe - 0.5) * (fo - 0.5) < 0:
                b = -b; action += "_eafswap"
        rows.append({**r.to_dict(), "beta_out_aligned": b, "harmonise": action})
    h = pd.DataFrame(rows)
    if len(h):
        h = h.rename(columns={"ea_exp": "ea", "nea_exp": "nea"})
        keep = ["rsid", "chr", "pos", "ea", "nea", "eaf", "info",
                "beta", "se", "p", "n", "F",
                "eaf_out", "beta_out_aligned", "se_out", "p_out", "harmonise", "dist_to_gene"]
        h = h[[c for c in keep if c in h.columns]]
    h.to_csv(f"{DIR}/mrinput_{out_name}.csv", index=False)
    print(f"  mrinput_{out_name}.csv: {len(h)} SNPs usable")
    return h

print("=== Step 1: clump ===")
instruments = {}
for ds, gene, core in CANDIDATES:
    inst = run_clump(ds, gene, core)
    if inst is not None and len(inst):
        instruments[(ds, gene, core)] = inst

print("\n=== Step 2: harmonize vs iPSYCH Model1/2 ===")
for (ds, gene, core), inst in instruments.items():
    tag = f"{ds}_{gene}_{'100kb' if core else '250kb'}"
    for model in ["Model1", "Model2", "Model3"]:
        opath = f"{DIR}/iPSYCH_{model}_{gene}_window.csv"
        if os.path.exists(opath):
            harmonize(inst, opath, f"{tag}_vs_iPSYCH_{model}")
print("\ndone")
