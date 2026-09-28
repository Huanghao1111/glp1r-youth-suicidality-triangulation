#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Child-side cis-MR: Wald ratio (1 SNP) / IVW (>=2 SNP) on harmonized mrinput tables.
Orientation: effect allele = BMI-LOWERING allele (drug-mimicking direction for GLP-1RA/MC4R agonism).
Outputs: child_mr_results.csv + child_mr_forest.png in E:/CM/GLP1/mr_input and workspace results/mr_cis."""
import glob, os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIRS = ["E:/CM/GLP1/mr_input", "results/mr_cis"]
files = sorted(glob.glob("E:/CM/GLP1/mr_input/mrinput_*.csv"))

def z2p(z):
    from math import erfc
    return erfc(abs(z) / np.sqrt(2))

rows, per_snp = [], []
for f in files:
    tag = os.path.basename(f).replace("mrinput_", "").replace(".csv", "")
    m = re.match(r"(.+)_(GLP1R|GIPR|MC4R)_(100kb|250kb)_vs_iPSYCH_(Model\d)", tag)
    if not m:
        continue
    ds, gene, win, model = m.groups()
    d = pd.read_csv(f)
    if len(d) == 0:
        continue
    # orient: effect allele = BMI-INCREASING allele -> theta = effect per +1 SD BMI
    flip = d["beta"] < 0
    d.loc[flip, "beta"] = -d.loc[flip, "beta"]
    d.loc[flip, "beta_out_aligned"] = -d.loc[flip, "beta_out_aligned"]
    # Wald ratio per SNP
    d["theta"] = d["beta_out_aligned"] / d["beta"]
    d["se_theta"] = d["se_out"] / d["beta"].abs()
    d["or_snp"] = np.exp(d["theta"])
    # Steiger (approx variance explained)
    d["r2_exp"] = 2 * d["eaf"] * (1 - d["eaf"]) * d["beta"] ** 2
    d["r2_out"] = 2 * d["eaf_out"] * (1 - d["eaf_out"]) * d["beta_out_aligned"] ** 2
    for _, r in d.iterrows():
        per_snp.append(dict(dataset=ds, gene=gene, window=win, model=model,
                            rsid=r["rsid"], beta_exp=r["beta"], theta=r["theta"],
                            or_snp=r["or_snp"], p_exp=r["p"], p_out=r["p_out"],
                            steiger_ok=r["r2_exp"] > r["r2_out"]))
    # IVW (fixed effect)
    w = (d["beta"] ** 2) / (d["se_out"] ** 2)
    theta = float((w * d["theta"]).sum() / w.sum())
    se = float(np.sqrt(1 / w.sum()))
    lo, hi = theta - 1.96 * se, theta + 1.96 * se
    # heterogeneity
    Q = float((w * (d["theta"] - theta) ** 2).sum())
    dfree = len(d) - 1
    p_het = z2p(np.sqrt(max(Q - 0, 0))) if False else (1 - 0)  # placeholder, compute properly below
    from scipy.stats import chi2
    p_het = 1 - chi2.cdf(Q, dfree) if dfree >= 1 else np.nan
    rows.append(dict(dataset=ds, gene=gene, window=win, model=model, n_snp=len(d),
                     snps=",".join(d["rsid"]), beta_MR=theta, se=se,
                     or_=np.exp(theta), or_lo=np.exp(lo), or_hi=np.exp(hi),
                     or_bmi_down=np.exp(-theta), or_bmi_down_lo=np.exp(-hi), or_bmi_down_hi=np.exp(-lo),
                     p=z2p(theta / se), Q=Q, p_het=p_het,
                     steiger_all_ok=bool((d["r2_exp"] > d["r2_out"]).all())))

res = pd.DataFrame(rows).sort_values(["gene", "dataset", "model"])
pd.set_option("display.width", 250)
show = res.copy()
show["OR (95%CI)"] = show.apply(lambda r: f"{r['or_']:.2f} ({r['or_lo']:.2f}-{r['or_hi']:.2f})", axis=1)
show["OR_BMI下降(药效向)"] = show.apply(lambda r: f"{r['or_bmi_down']:.2f} ({r['or_bmi_down_lo']:.2f}-{r['or_bmi_down_hi']:.2f})", axis=1)
print(show[["dataset", "gene", "model", "n_snp", "OR (95%CI)", "OR_BMI下降(药效向)", "p", "p_het", "steiger_all_ok"]].to_string(index=False))
res.to_csv("E:/CM/GLP1/mr_input/child_mr_results.csv", index=False)
res.to_csv("results/mr_cis/child_mr_results.csv", index=False)
pd.DataFrame(per_snp).to_csv("E:/CM/GLP1/mr_input/child_mr_per_snp.csv", index=False)
pd.DataFrame(per_snp).to_csv("results/mr_cis/child_mr_per_snp.csv", index=False)

# ---- forest plot (GLP1R only; MC4R/GIPR dropped from main line) ----
plt.rcParams["font.family"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
fig, axes = plt.subplots(1, 1, figsize=(8.5, 6.5))
axes = [axes]
for ax, gene in zip(axes, ["GLP1R"]):
    sub = res[res["gene"] == gene].copy()
    if len(sub) == 0:
        continue
    sub["label"] = [
        f"{d.replace('MoBa2022_','').replace('EGG2020_childBMI','EGG(2-10岁)')} | {md.replace('Model','M')} (n={n})"
        for d, md, n in zip(sub["dataset"], sub["model"], sub["n_snp"])
    ]
    sub = sub.iloc[::-1]
    y = np.arange(len(sub))
    ax.errorbar(sub["or_"], y, xerr=[sub["or_"] - sub["or_lo"], sub["or_hi"] - sub["or_"]],
                fmt="o", color="#2166ac", ecolor="#2166ac", elinewidth=1.6, capsize=3, ms=5)
    ax.axvline(1, color="grey", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_yticks(y); ax.set_yticklabels(sub["label"], fontsize=8.5)
    ax.set_xlabel("OR per +1 SD childhood BMI (log scale)；OR<1 = BMI越高风险越低")
    ax.set_title(f"{gene} 顺式工具 → iPSYCH 自杀未遂", fontsize=11)
    for yi, (_, r) in zip(y, sub.iterrows()):
        ax.text(ax.get_xlim()[1] * 1.1, yi, f"{r['or_']:.2f} ({r['or_lo']:.2f}-{r['or_hi']:.2f})",
                va="center", fontsize=8)
    ax.set_xlim(right=ax.get_xlim()[1] * 3.2)
fig.suptitle("儿童期 BMI（GLP1R 顺式工具）→ 青少年自杀未遂（iPSYCH，<32 岁）\n注意：药效方向 = BMI 下降 = 取倒数（OR>1 = 风险升高）", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
for out in ["E:/CM/GLP1/mr_input/child_mr_forest.png", "results/mr_cis/child_mr_forest.png"]:
    fig.savefig(out, dpi=200, bbox_inches="tight")
print("\nsaved child_mr_results.csv / child_mr_per_snp.csv / child_mr_forest.png")
