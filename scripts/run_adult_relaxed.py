#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Adult GLP1R relaxed-threshold scan (ieu-b-40, +-100kb window).
Main-line p<5e-8 yields 0 candidates; here we relax p in a ladder,
clump with the SAME parameters as the main line (PLINK r2=0.3, kb=100, 1000G EUR),
harmonize vs iPSYCH Model1/2/3 and run fixed-effect IVW (Wald if 1 SNP).
Exploratory sensitivity analysis: sub-GWS instruments, winner's curse caveat.
Output: adult_glp1r_relaxed.csv + mrinput tables per threshold (E:/CM/GLP1/mr_input,
workspace results/mr_cis backup)."""
import os, subprocess
import numpy as np
import pandas as pd
from scipy.stats import chi2
from math import erfc

DIR = "E:/CM/GLP1/mr_input"
WS = "results/mr_cis"
PLINK = "E:/CM/GLP1/LD/plink_win64/plink.exe"
BFILE = "E:/CM/GLP1/LD/bfile/EUR"
PALINDROMIC = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
THRESHOLDS = [5e-7, 5e-6, 1e-5, 1e-4]

def z2p(z):
    return erfc(abs(z) / np.sqrt(2))

def run_clump(cand, thr, tag):
    inp = f"{DIR}/clumpin_{tag}.txt"
    cand[["rsid", "p"]].rename(columns={"rsid": "SNP", "p": "P"}).to_csv(inp, sep=" ", index=False)
    outpref = f"{DIR}/clump_{tag}"
    cmd = [PLINK, "--bfile", BFILE, "--clump", inp, "--clump-p1", str(thr),
           "--clump-p2", "1", "--clump-r2", "0.3", "--clump-kb", "100",
           "--clump-snp-field", "SNP", "--clump-field", "P", "--out", outpref, "--silent"]
    subprocess.run(cmd, capture_output=True, text=True)
    cl = f"{outpref}.clumped"
    if not os.path.exists(cl):
        return []
    tbl = pd.read_csv(cl, sep=r"\s+", engine="python")
    return tbl["SNP"].tolist()

def harmonize(inst, outcome_path):
    out = pd.read_csv(outcome_path)
    o = out.rename(columns={"beta": "beta_out", "se": "se_out", "p": "p_out", "eaf": "eaf_out"})
    o = o[["rsid", "chr", "pos", "ea", "nea", "eaf_out", "beta_out", "se_out", "p_out"]]
    m = inst.merge(o.drop(columns=["chr", "pos"]), on="rsid", how="inner", suffixes=("_exp", "_outc"))
    rest = inst[~inst["rsid"].isin(set(m["rsid"]))]
    if len(rest):
        o2 = o.rename(columns={"rsid": "rsid_out"})
        m2 = rest.merge(o2, on=["chr", "pos"], how="inner", suffixes=("_exp", "_outc"))
        if len(m2):
            m2["rsid"] = m2["rsid"] + "|" + m2["rsid_out"].astype(str)
            m2 = m2.drop(columns=["rsid_out"])
            m = pd.concat([m, m2[m.columns]], ignore_index=True)
    rows = []
    for _, r in m.iterrows():
        ea_e, nea_e = str(r["ea_exp"]).upper(), str(r["nea_exp"]).upper()
        ea_o, nea_o = str(r["ea_outc"]).upper(), str(r["nea_outc"]).upper()
        if (ea_o, nea_o) == (ea_e, nea_e):
            b = r["beta_out"]; action = "aligned"
        elif (ea_o, nea_o) == (nea_e, ea_e):
            b = -r["beta_out"]; action = "flipped"
        else:
            continue
        if (ea_e, nea_e) in PALINDROMIC:
            fe, fo = r.get("eaf_exp", np.nan), r.get("eaf_out", np.nan)
            if pd.isna(fe) or pd.isna(fo) or (0.4 < fe < 0.6) or (0.4 < fo < 0.6):
                continue
            if (fe - 0.5) * (fo - 0.5) < 0:
                b = -b; action += "_eafswap"
        rows.append({**r.to_dict(), "beta_out_aligned": b, "harmonise": action})
    return pd.DataFrame(rows)

def ivw(h):
    flip = h["beta"] < 0
    h = h.copy()
    h.loc[flip, "beta"] = -h.loc[flip, "beta"]
    h.loc[flip, "beta_out_aligned"] = -h.loc[flip, "beta_out_aligned"]
    h["theta"] = h["beta_out_aligned"] / h["beta"]
    w = (h["beta"] ** 2) / (h["se_out"] ** 2)
    theta = float((w * h["theta"]).sum() / w.sum())
    se = float(np.sqrt(1 / w.sum()))
    Q = float((w * (h["theta"] - theta) ** 2).sum()) if len(h) > 1 else np.nan
    return dict(beta_MR=theta, se=se, or_=np.exp(theta),
                or_lo=np.exp(theta - 1.96 * se), or_hi=np.exp(theta + 1.96 * se),
                p=z2p(theta / se),
                p_het=(1 - chi2.cdf(Q, len(h) - 1)) if len(h) > 1 else np.nan), h

win = pd.read_csv(f"{DIR}/AdultBMI_ieub40_GLP1R_window.csv")
core = win[win["in_100kb"]].copy()
core = core[core["rsid"].astype(str).str.startswith("rs")]

all_rows = []
for thr in THRESHOLDS:
    cand = core[core["p"] < thr].copy()
    tag = f"AdultBMI_ieub40_GLP1R_100kb_p{thr:.0e}".replace("-", "m")
    if len(cand) == 0:
        print(f"p<{thr:.0e}: 0 candidates"); continue
    idx = run_clump(cand, thr, tag)
    inst = cand[cand["rsid"].isin(idx)].copy()
    inst.to_csv(f"{DIR}/instruments_{tag}.csv", index=False)
    print(f"\np<{thr:.0e}: {len(cand)} candidates -> {len(idx)} independent: {idx}")
    if len(idx) == 0:
        continue
    meanF = float(((inst["beta"] / inst["se"]) ** 2).mean())
    minF = float(((inst["beta"] / inst["se"]) ** 2).min())
    for model in ["Model1", "Model2", "Model3"]:
        opath = f"{DIR}/iPSYCH_{model}_GLP1R_window.csv"
        h = harmonize(inst, opath)
        h.to_csv(f"{DIR}/mrinput_{tag}_vs_iPSYCH_{model}.csv", index=False)
        if len(h) == 0:
            print(f"  {model}: 0 harmonized"); continue
        r, hdetail = ivw(h)
        all_rows.append(dict(threshold=thr, model=model, n_cand=len(cand), n_indep=len(idx),
                             n_harm=len(h), mean_F=meanF, min_F=minF, **r))
        print(f"  {model}: n={len(h)}  OR(+1SD)={r['or_']:.3f} ({r['or_lo']:.3f}-{r['or_hi']:.3f}) "
              f"p={r['p']:.4g}  OR(-1SD药效)={1/r['or_']:.3f}  meanF={meanF:.1f}")

out = pd.DataFrame(all_rows)
out["or_down"] = 1 / out["or_"]
out["or_down_lo"] = 1 / out["or_hi"]
out["or_down_hi"] = 1 / out["or_lo"]
out.to_csv(f"{DIR}/adult_glp1r_relaxed.csv", index=False)
out.to_csv(f"{WS}/adult_glp1r_relaxed.csv", index=False)
print("\nsaved adult_glp1r_relaxed.csv")
