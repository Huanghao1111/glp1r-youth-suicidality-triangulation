#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Steiger-filtered cis-MR sensitivity (standard TwoSampleMR steiger_filtering).
Per-SNP Steiger directionality with formal Fisher-z test (independent samples):
  r_exp = beta * sqrt(2*eaf*(1-eaf))            (BMI, standardized)
  r_out = z/sqrt(n_eff-1+z^2), z=beta/se        (iPSYCH case-control)
  z_f   = (atanh(r_exp)-atanh(r_out))/sqrt(1/(n_exp-3)+1/(n_out-3))
  drop  = r_out^2 > r_exp^2 AND P < .05         (significantly WRONG direction only;
          correct-but-unsignificant variants are kept, per TwoSampleMR default)
IVW (fixed effect, same as run_child_mr.py) re-run on kept variants.
Scope: GLP1R, five MoBa infant waves + adult BMI meta, iPSYCH Models 1-3.
Output: results/mr_cis/steiger_filter_sensitivity.csv (+ _per_snp.csv)"""
import glob, os, re
import numpy as np
import pandas as pd
from math import erfc, atanh

IPSYCH_NEFF = {"Model1": 4*6024*44240/(6024+44240),
               "Model2": 4*6024*44240/(6024+44240),
               "Model3": 4*4302*13294/(4302+13294)}

def z2p(z):
    return erfc(abs(z) / np.sqrt(2))

def ivw(d):
    w = d["beta"] ** 2 / d["se_out"] ** 2
    th = float((w * d["theta"]).sum() / w.sum())
    se = float(np.sqrt(1 / w.sum()))
    return th, se, z2p(th / se)

rows, per = [], []
pats = []
for f in sorted(glob.glob("E:/CM/GLP1/mr_input/mrinput_*_GLP1R_100kb_vs_iPSYCH_*.csv")):
    tag = os.path.basename(f)[8:-4]
    m = re.match(r"(.+)_GLP1R_100kb_vs_iPSYCH_(Model\d)", tag)
    ds, model = m.group(1), m.group(2)
    if not (ds.startswith("MoBa2022") or ds == "AdultBMI_meta"):
        continue
    d = pd.read_csv(f)
    if not len(d):
        continue
    flip = d["beta"] < 0
    d.loc[flip, "beta"] = -d.loc[flip, "beta"]
    d.loc[flip, "beta_out_aligned"] = -d.loc[flip, "beta_out_aligned"]
    d["theta"] = d["beta_out_aligned"] / d["beta"]
    d["se_theta"] = d["se_out"] / d["beta"].abs()
    n_out = IPSYCH_NEFF[model]
    n_exp = d["n"] if "n" in d and d["n"].notna().all() else pd.Series(np.nanmedian(d["n"]), index=d.index)
    r_exp = d["beta"] * np.sqrt(2 * d["eaf"] * (1 - d["eaf"]))
    z_o = d["beta_out_aligned"] / d["se_out"]
    r_out = z_o / np.sqrt(n_out - 1 + z_o ** 2)
    zf = (np.arctanh(np.clip(r_exp, -0.9999, 0.9999)) - np.arctanh(np.clip(r_out, -0.9999, 0.9999))) / \
         np.sqrt(1 / (n_exp - 3) + 1 / (n_out - 3))
    d["r_exp"], d["r_out"] = r_exp, r_out
    d["steiger_z"], d["steiger_p"] = zf, [z2p(z) for z in zf]
    # standard TwoSampleMR steiger_filtering: drop only significantly WRONG direction
    d["steiger_pass"] = ~((d["r_out"] ** 2 > d["r_exp"] ** 2) & (d["steiger_p"] < 0.05))
    th0, se0, p0 = ivw(d)
    keep = d[d["steiger_pass"]]
    if len(keep) >= 1:
        th1, se1, p1 = ivw(keep)
        or1, lo1, hi1 = np.exp(th1), np.exp(th1 - 1.96 * se1), np.exp(th1 + 1.96 * se1)
    else:
        or1 = lo1 = hi1 = p1 = np.nan
    rows.append(dict(dataset=ds, model=model, n_snp_all=len(d), n_snp_kept=len(keep),
                     dropped=",".join(d.loc[~d["steiger_pass"], "rsid"]) or "—",
                     or_all=np.exp(th0), or_all_lo=np.exp(th0 - 1.96 * se0), or_all_hi=np.exp(th0 + 1.96 * se0), p_all=p0,
                     or_steiger=or1, or_lo=lo1, or_hi=hi1, p_steiger=p1))
    for _, r in d.iterrows():
        per.append(dict(dataset=ds, model=model, rsid=r["rsid"], r_exp=round(r["r_exp"], 5),
                        r_out=round(r["r_out"], 5), steiger_z=round(r["steiger_z"], 2),
                        steiger_p=float(r["steiger_p"]), passed=bool(r["steiger_pass"])))

res = pd.DataFrame(rows)
pd.set_option("display.width", 250)
show = res.copy()
for c in ["or_all", "or_all_lo", "or_all_hi", "or_steiger", "or_lo", "or_hi"]:
    show[c] = show[c].astype(float).round(3)
show["p_all"] = show["p_all"].map(lambda x: f"{x:.2g}")
show["p_steiger"] = show["p_steiger"].map(lambda x: f"{x:.2g}" if pd.notna(x) else "—")
print(show.to_string(index=False))
res.to_csv("results/mr_cis/steiger_filter_sensitivity.csv", index=False)
res.to_csv("E:/CM/GLP1/mr_input/steiger_filter_sensitivity.csv", index=False)
pd.DataFrame(per).to_csv("results/mr_cis/steiger_filter_per_snp.csv", index=False)
print("\nsaved steiger_filter_sensitivity.csv + steiger_filter_per_snp.csv")
