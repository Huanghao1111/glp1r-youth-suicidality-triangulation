#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""coloc.abf (Giambartolomei 2014) implemented in Python, GLP1R +-250kb region.
Wakefield ABF per SNP (sign-agnostic; beta/se required for both traits):
  logABF = 0.5*log(V/(V+W^2)) + 0.5*(W^2/(V+W^2))*z^2,  V=se^2
Posteriors: H0=1; H1=p1*S1; H2=p2*S2; H3=p1*p2*(S1*S2 - S12); H4=p12*S12
with S1=sum ABF1, S2=sum ABF2, S12=sum ABF1*ABF2 (log-space, logsumexp).
Priors p1=p2=1e-4, p12=1e-5 (coloc defaults). W=0.15 quantitative (BMI),
W=0.20 case-control (iPSYCH). Single-causal-variant assumption caveat:
region harbours 2-3 near-independent BMI signals -> H4 may be deflated.
Output: coloc_glp1r.csv (E:/CM/GLP1/mr_input + results/mr_cis)."""
import numpy as np
import pandas as pd
from scipy.special import logsumexp

DIR = "E:/CM/GLP1/mr_input"
WS = "results/mr_cis"
P1, P2, P12 = 1e-4, 1e-4, 1e-5
W_Q, W_CC = 0.15, 0.20

def logabf(beta, se, W):
    V = se ** 2
    r = W ** 2 / (V + W ** 2)
    return 0.5 * np.log(V / (V + W ** 2)) + 0.5 * r * (beta / se) ** 2

def load_win(path):
    df = pd.read_csv(path)
    df = df.dropna(subset=["beta", "se"])
    df = df[(df["se"] > 0)]
    df = df.sort_values("p").drop_duplicates("rsid", keep="first")
    return df[["rsid", "chr", "pos", "beta", "se", "p"]]

def merge_pair(a, b):
    m = a.merge(b, on="rsid", how="inner", suffixes=("_1", "_2"))
    hit = set(m["rsid"])
    rest = a[~a["rsid"].isin(hit)]
    if len(rest):  # position fallback (exm probe ids in iPSYCH)
        b2 = b.rename(columns={"rsid": "rsid_2b"})
        m2 = rest.merge(b2, on=["chr", "pos"], how="inner", suffixes=("_1", "_2"))
        if len(m2):
            m2 = m2.assign(rsid=m2["rsid"] + "|" + m2["rsid_2b"].astype(str),
                           chr_1=m2["chr"], pos_1=m2["pos"],
                           chr_2=m2["chr"], pos_2=m2["pos"])
            m = pd.concat([m, m2[m.columns]], ignore_index=True)
    m = m.sort_values("pos_1").drop_duplicates(subset=["chr_1", "pos_1"], keep="first")
    return m

def coloc_abf(m):
    l1 = logabf(m["beta_1"].values, m["se_1"].values, W_Q)
    l2 = logabf(m["beta_2"].values, m["se_2"].values, W_CC)
    S1 = logsumexp(l1); S2 = logsumexp(l2); S12 = logsumexp(l1 + l2)
    # H3 needs S1*S2 - S12 in linear space: log-diff
    a = S1 + S2
    ld = a + np.log1p(-np.exp(S12 - a)) if S12 < a else -np.inf
    lw = np.array([0.0, np.log(P1) + S1, np.log(P2) + S2,
                   np.log(P1) + np.log(P2) + ld, np.log(P12) + S12])
    pp = np.exp(lw - logsumexp(lw))
    top = int(np.argmax(l1 + l2))
    return pp, m.iloc[top]

PAIRS = [
    ("MoBa2022_3months", "Model1"), ("MoBa2022_3months", "Model2"), ("MoBa2022_3months", "Model3"),
    ("MoBa2022_1year", "Model1"),
    ("AdultBMI_meta", "Model1"), ("AdultBMI_meta", "Model3"),
]
rows = []
for exp, model in PAIRS:
    a = load_win(f"{DIR}/{exp}_GLP1R_window.csv")
    b = load_win(f"{DIR}/iPSYCH_{model}_GLP1R_window.csv")
    m = merge_pair(a, b)
    pp, top = coloc_abf(m)
    rows.append(dict(exposure=exp, outcome=f"iPSYCH_{model}", n_snp=len(m),
                     PP_H0=pp[0], PP_H1=pp[1], PP_H2=pp[2], PP_H3=pp[3], PP_H4=pp[4],
                     top_shared_rsid=top["rsid"],
                     top_p_exp=top["p_1"], top_p_out=top["p_2"]))
    print(f"{exp:20s} x {model}: n={len(m):4d}  H3={pp[3]:.3f}  H4={pp[4]:.3f}  "
          f"top={top['rsid']} (p_exp={top['p_1']:.1e}, p_out={top['p_2']:.1e})")

out = pd.DataFrame(rows)
out.to_csv(f"{DIR}/coloc_glp1r.csv", index=False)
out.to_csv(f"{WS}/coloc_glp1r.csv", index=False)
print("\nsaved coloc_glp1r.csv")
