#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""coloc.susie-style analysis, GLP1R +-250kb region.
SuSiE-RSS (IBSS with summary stats + external LD) implemented in numpy:
  - z-scores harmonised to 1000G EUR bim A1 allele (sign flip when ea==A2)
  - LD matrix from EUR panel via PLINK --r square (per exposure-outcome pair,
    on the pair's merged SNP set, same SNPs both traits -> coloc.susie logic)
  - per-effect prior variance estimated by grid on SER marginal likelihood
  - 95% credible sets; purity filter (min |r| >= 0.5), near-uniform effects dropped
Then coloc.abf (same Wakefield ABF + priors as run_coloc.py) for every pair of
credible sets (intersection >= 5 SNPs), reporting the best-H4 pair per
exposure x outcome, alongside the previous single-variant coloc.abf result.
Output: results/mr_cis/coloc_susie_glp1r_summary.csv
        results/mr_cis/coloc_susie_glp1r_pairs.csv
        results/mr_cis/coloc_susie_glp1r_cs.csv"""
import os, subprocess, tempfile
import numpy as np
import pandas as pd
from scipy.special import logsumexp

DIR = "E:/CM/GLP1/mr_input"
BFILE = "E:/CM/GLP1/LD/bfile/EUR"
PLINK = "E:/CM/GLP1/LD/plink_win64/plink.exe"
WS = "results/mr_cis"
P1, P2, P12 = 1e-4, 1e-4, 1e-5
W_Q, W_CC = 0.15, 0.20
VGRID = np.array([1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1])
TMP = tempfile.mkdtemp(prefix="colocsusie_")

# ---------- data loading / harmonisation ----------
def load_bim():
    bim = pd.read_csv(BFILE + ".bim", sep="\t", header=None,
                      names=["chr", "rsid", "cm", "pos", "A1", "A2"])
    return bim.set_index("rsid")

BIM = load_bim()

def load_win(path, n_fix=None):
    df = pd.read_csv(path)
    df = df.dropna(subset=["beta", "se"])
    df = df[df["se"] > 0].sort_values("p").drop_duplicates("rsid", keep="first")
    if "n" not in df or df["n"].isna().all():
        df["n"] = np.nan
    df["n"] = df["n"].fillna(n_fix if n_fix else np.nan)
    df = df.dropna(subset=["n"])
    return df[["rsid", "chr", "pos", "ea", "nea", "beta", "se", "p", "n"]]

def harmonise(df):
    """flip beta to bim A1 orientation; drop SNPs not in panel or allele-mismatched"""
    d = df.merge(BIM[["A1", "A2"]], left_on="rsid", right_index=True, how="inner")
    ok = ((d["ea"] == d["A1"]) & (d["nea"] == d["A2"])) | \
         ((d["ea"] == d["A2"]) & (d["nea"] == d["A1"]))
    d = d[ok].copy()
    flip = d["ea"] == d["A2"]
    d.loc[flip, "beta"] = -d.loc[flip, "beta"]
    d["z"] = d["beta"] / d["se"]
    return d

def merge_pair(a, b):
    m = a.merge(b, on="rsid", how="inner", suffixes=("_1", "_2"))
    return m.sort_values("pos_1").reset_index(drop=True)

# ---------- LD via PLINK ----------
def ld_matrix(rsids, tag):
    lst = os.path.join(TMP, f"{tag}.snps")
    with open(lst, "w") as f:
        f.write("\n".join(rsids))
    out = os.path.join(TMP, tag)
    subprocess.run([PLINK, "--bfile", BFILE, "--extract", lst,
                    "--r", "square", "--write-snplist", "--out", out,
                    "--silent"], check=True)
    order = [l.strip() for l in open(out + ".snplist")]
    R = np.loadtxt(out + ".ld")
    if R.ndim == 1:  # single SNP edge
        R = R.reshape(1, 1)
    return order, R

# ---------- SuSiE-RSS ----------
def ser(r, V, S):
    """single-effect regression; r = residualised Xty on the (n-1) scale,
    XtX_jj = S. returns alpha, mu, marginal log-evidence"""
    s2 = V / (1 + S * V)
    lbf = -0.5 * np.log1p(S * V) + 0.5 * V * r ** 2 / (1 + S * V)
    la = np.log(1 / len(r)) + lbf
    la -= la.max()
    a = np.exp(la); a /= a.sum()
    return a, s2 * r, logsumexp(np.log(1 / len(r)) + lbf)

def susie_rss(z, n, R, L=10, max_iter=300, tol=1e-7):
    """z + per-trait median n -> correlation scale, then Xty = S*cor, XtX = S*R"""
    S = float(np.median(n)) - 1.0
    Xty = S * z / np.sqrt(S + z ** 2)
    XtX = S * R
    p = len(z)
    alpha = np.full((L, p), 1.0 / p)
    mu = np.zeros((L, p))
    Vs = np.full(L, 0.001)
    total = np.zeros(p)
    for it in range(max_iter):
        prev = total.copy()
        for l in range(L):
            total = total - alpha[l] * mu[l]
            r = Xty - XtX @ total
            best = None
            for V in VGRID:
                a, m_, ev = ser(r, V, S)
                if best is None or ev > best[2]:
                    best = (a, m_, ev, V)
            alpha[l], mu[l], Vs[l] = best[0], best[1], best[3]
            total = total + alpha[l] * mu[l]
        if np.max(np.abs(total - prev)) < tol:
            break
    # prune effects whose final SER evidence is worse than null
    keep_l = []
    for l in range(L):
        total_minus = total - alpha[l] * mu[l]
        r = Xty - XtX @ total_minus
        ev = ser(r, Vs[l], S)[2]
        if ev > 0:
            keep_l.append(l)
    pip = 1 - np.prod(1 - alpha[keep_l], axis=0) if keep_l else np.zeros(p)
    cs = []
    for l in keep_l:
        order = np.argsort(-alpha[l])
        cum = np.cumsum(alpha[l][order])
        k = int(np.searchsorted(cum, 0.95)) + 1
        idx = np.sort(order[:k])
        purity = 1.0 if len(idx) == 1 else np.min(np.abs(R[np.ix_(idx, idx)]))
        cs.append(dict(effect=l, idx=idx, size=len(idx), purity=purity,
                       max_alpha=float(alpha[l].max()), V=Vs[l]))
    return dict(alpha=alpha, mu=mu, pip=pip, cs=cs, Vs=Vs, iters=it + 1)

# ---------- coloc.abf (unchanged from run_coloc.py) ----------
def logabf(beta, se, W):
    V = se ** 2
    r = W ** 2 / (V + W ** 2)
    return 0.5 * np.log(V / (V + W ** 2)) + 0.5 * r * (beta / se) ** 2

def coloc_abf(b1, s1, b2, s2):
    l1 = logabf(b1, s1, W_Q); l2 = logabf(b2, s2, W_CC)
    S1 = logsumexp(l1); S2 = logsumexp(l2); S12 = logsumexp(l1 + l2)
    a = S1 + S2
    ld = a + np.log1p(-np.exp(S12 - a)) if S12 < a else -np.inf
    lw = np.array([0.0, np.log(P1) + S1, np.log(P2) + S2,
                   np.log(P1) + np.log(P2) + ld, np.log(P12) + S12])
    return np.exp(lw - logsumexp(lw))

# ---------- main ----------
# iPSYCH effective n = 4/(1/ncase + 1/nctrl); S15 footnote:
# Model 1/2: 6,024 cases, 44,240 controls; Model 3: 4,302 cases, 13,294 controls
def n_eff(nc, nt):
    return 4 * nc * nt / (nc + nt)

IPSYCH_N = {"Model1": n_eff(6024, 44240), "Model2": n_eff(6024, 44240),
            "Model3": n_eff(4302, 13294)}

PAIRS = [
    ("MoBa2022_3months", "Model1"), ("MoBa2022_3months", "Model2"), ("MoBa2022_3months", "Model3"),
    ("MoBa2022_1year", "Model1"),
    ("AdultBMI_meta", "Model1"), ("AdultBMI_meta", "Model3"),
]
summary, pair_rows, cs_rows = [], [], []
for exp, model in PAIRS:
    tag = f"{exp}_{model}"
    a = harmonise(load_win(f"{DIR}/{exp}_GLP1R_window.csv"))
    b = harmonise(load_win(f"{DIR}/iPSYCH_{model}_GLP1R_window.csv",
                           n_fix=IPSYCH_N[model]))
    m = merge_pair(a, b)
    order, R = ld_matrix(list(m["rsid"]), tag)
    keep = m["rsid"].isin(order)
    m = m[keep].reset_index(drop=True)
    idx = [order.index(r) for r in m["rsid"]]
    R = R[np.ix_(idx, idx)]
    np.fill_diagonal(R, 1.0)
    R = np.clip(R, -0.999, 0.999)

    fits = []
    for sfx in ("_1", "_2"):
        fits.append(susie_rss(m[f"z{sfx}"].values, m[f"n{sfx}"].values, R))
    print(f"\n=== {exp} x iPSYCH_{model}: n={len(m)} SNPs ===")
    for t, fit in enumerate(fits):
        good = [c for c in fit["cs"] if c["purity"] >= 0.5 and c["size"] <= 0.95 * len(m)
                and c["max_alpha"] > 3.0 / len(m)]
        fits[t]["good_cs"] = good
        print(f"  trait{t+1}: {len(good)} usable CS (of {len(fit['cs'])} effects), "
              f"maxPIP={fit['pip'].max():.3f}, iters={fit['iters']}")
        for c in good:
            cs_rows.append(dict(pair=tag, trait=t + 1, effect=c["effect"],
                                size=c["size"], purity=round(c["purity"], 3),
                                max_alpha=round(c["max_alpha"], 4),
                                top_rsid=m["rsid"].iloc[c["idx"][np.argmax(fits[t]['alpha'][c['effect']][c['idx']])]]))
    # coloc over CS pairs
    best = None
    for c1 in fits[0].get("good_cs", []):
        for c2 in fits[1].get("good_cs", []):
            inter = np.intersect1d(c1["idx"], c2["idx"])
            if len(inter) < 5:
                continue
            pp = coloc_abf(m["beta_1"].values[inter], m["se_1"].values[inter],
                           m["beta_2"].values[inter], m["se_2"].values[inter])
            row = dict(pair=tag, cs_exp=c1["effect"], cs_out=c2["effect"],
                       n_shared=len(inter), PP_H3=pp[3], PP_H4=pp[4],
                       top_rsid=m["rsid"].iloc[inter[np.argmax(
                           logabf(m["beta_1"].values[inter], m["se_1"].values[inter], W_Q) +
                           logabf(m["beta_2"].values[inter], m["se_2"].values[inter], W_CC))]])
            pair_rows.append(row)
            if best is None or pp[4] > best["PP_H4"]:
                best = row
            print(f"    CS{c1['effect']}xCS{c2['effect']} n={len(inter):3d} "
                  f"H3={pp[3]:.3f} H4={pp[4]:.3f}")
    summary.append(dict(exposure=exp, outcome=f"iPSYCH_{model}", n_snp=len(m),
                        n_cs_exp=len(fits[0].get("good_cs", [])),
                        n_cs_out=len(fits[1].get("good_cs", [])),
                        best_PP_H3=best["PP_H3"] if best else np.nan,
                        best_PP_H4=best["PP_H4"] if best else np.nan,
                        best_pair=f"CS{best['cs_exp']}xCS{best['cs_out']}" if best else "—"))

s = pd.DataFrame(summary)
s.to_csv(f"{WS}/coloc_susie_glp1r_summary.csv", index=False)
pd.DataFrame(pair_rows).to_csv(f"{WS}/coloc_susie_glp1r_pairs.csv", index=False)
pd.DataFrame(cs_rows).to_csv(f"{WS}/coloc_susie_glp1r_cs.csv", index=False)
print("\n=== summary (best CS pair per exposure-outcome) ===")
print(s.to_string(index=False))
print("\nsaved coloc_susie_glp1r_{summary,pairs,cs}.csv; tmp LD in", TMP)
