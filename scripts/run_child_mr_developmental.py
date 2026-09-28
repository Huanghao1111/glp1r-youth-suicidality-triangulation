#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Fixed-instrument developmental scan, both genes:
GLP1R union instruments (infancy-selected) + MC4R union instruments (rs11873305 MoBa-7y,
rs953442 EGG downstream) looked up across ALL 12 MoBa timepoints -> IVW vs iPSYCH M1/2/3.
Orientation: per +1 SD BMI. Winner's curse caveat: instruments selected on their
respective significant timepoint."""
import glob, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2
from math import erfc

DIR = "E:/CM/GLP1/mr_input"
WS = "results/mr_cis"
TPS = ["birth", "6weeks", "3months", "6months", "8months", "1year", "1.5years",
       "2years", "3years", "5years", "7years", "8years"]
AGE = {"birth": 0, "6weeks": 0.12, "3months": 0.25, "6months": 0.5, "8months": 0.67,
       "1year": 1, "1.5years": 1.5, "2years": 2, "3years": 3, "5years": 5,
       "7years": 7, "8years": 8}
PAL = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

def z2p(z): return erfc(abs(z) / np.sqrt(2))

def harmonize(exp, outc):
    o = outc.rename(columns={"beta": "beta_out", "se": "se_out", "p": "p_out", "eaf": "eaf_out"})
    m = exp.merge(o[["rsid", "ea", "nea", "eaf_out", "beta_out", "se_out", "p_out"]],
                  on="rsid", how="inner", suffixes=("_exp", "_outc"))
    hit = set(m["rsid"])
    rest = exp[~exp["rsid"].isin(hit)]
    if len(rest):  # position fallback (exm ids in iPSYCH)
        o2 = outc.rename(columns={"beta": "beta_out", "se": "se_out", "p": "p_out", "eaf": "eaf_out",
                                  "rsid": "rsid_out"})
        m2 = rest.merge(o2[["rsid_out", "chr", "pos", "ea", "nea", "eaf_out", "beta_out", "se_out", "p_out"]],
                        on=["chr", "pos"], how="inner", suffixes=("_exp", "_outc"))
        if len(m2):
            m2["rsid"] = m2["rsid"].astype(str) + "|" + m2["rsid_out"].astype(str)
            m = pd.concat([m, m2[m.columns]], ignore_index=True)
    rows = []
    for _, r in m.iterrows():
        ea_e, nea_e = str(r["ea_exp"]).upper(), str(r["nea_exp"]).upper()
        ea_o, nea_o = str(r["ea_outc"]).upper(), str(r["nea_outc"]).upper()
        if (ea_o, nea_o) == (ea_e, nea_e): b = r["beta_out"]
        elif (ea_o, nea_o) == (nea_e, ea_e): b = -r["beta_out"]
        else: continue
        if (ea_e, nea_e) in PAL:
            fe, fo = r.get("eaf_exp", np.nan), r.get("eaf_out", np.nan)
            if pd.isna(fe) or pd.isna(fo) or (0.4 < fe < 0.6) or (0.4 < fo < 0.6): continue
            if (fe - 0.5) * (fo - 0.5) < 0: b = -b
        rows.append({**r.to_dict(), "beta_out_aligned": b})
    return pd.DataFrame(rows)

all_res = []
for gene, inst_glob in [("GLP1R", "instruments_MoBa2022_*_GLP1R_100kb.csv")]:
    if gene == "MC4R":
        inst = pd.concat([pd.read_csv(f) for f in
                          glob.glob(f"{DIR}/instruments_MoBa2022_*_MC4R_100kb.csv") +
                          glob.glob(f"{DIR}/instruments_EGG2020_childBMI_MC4R_250kb.csv")])
    else:
        inst = pd.concat([pd.read_csv(f) for f in glob.glob(f"{DIR}/{inst_glob}")])
    union = inst.drop_duplicates("rsid")[["rsid", "chr", "pos", "ea", "nea"]]
    print(f"{gene} union instruments: {len(union)} -> {union['rsid'].tolist()}")
    outcomes = {m: pd.read_csv(f"{DIR}/iPSYCH_{m}_{gene}_window.csv") for m in ["Model1", "Model2", "Model3"]}
    for tp in TPS:
        win = pd.read_csv(f"{DIR}/MoBa2022_{tp}_{gene}_window.csv")
        exp = win[win["rsid"].isin(union["rsid"])].copy()
        if len(exp) == 0: continue
        meanF = float(((exp["beta"] / exp["se"]) ** 2).mean())
        minp = float(win[win["in_100kb"]]["p"].min())
        for model, outc in outcomes.items():
            h = harmonize(exp, outc)
            if len(h) < 2:
                all_res.append(dict(gene=gene, timepoint=tp, age=AGE[tp], model=model,
                                    n_snp=len(h), mean_F=meanF, min_p_100kb=minp)); continue
            flip = h["beta"] < 0   # orient: per +1 SD BMI
            h.loc[flip, "beta"] = -h.loc[flip, "beta"]
            h.loc[flip, "beta_out_aligned"] = -h.loc[flip, "beta_out_aligned"]
            h["theta"] = h["beta_out_aligned"] / h["beta"]
            w = (h["beta"] ** 2) / (h["se_out"] ** 2)
            theta = float((w * h["theta"]).sum() / w.sum())
            se = float(np.sqrt(1 / w.sum()))
            Q = float((w * (h["theta"] - theta) ** 2).sum())
            all_res.append(dict(gene=gene, timepoint=tp, age=AGE[tp], model=model, n_snp=len(h),
                                beta_MR=theta, se=se, or_=np.exp(theta),
                                or_lo=np.exp(theta - 1.96 * se), or_hi=np.exp(theta + 1.96 * se),
                                p=z2p(theta / se), p_het=1 - chi2.cdf(Q, len(h) - 1),
                                mean_F=meanF, min_p_100kb=minp))

res = pd.DataFrame(all_res)

# ---- adult pseudo-timepoint: GLP1R fixed infancy union; MC4R own adult instruments ----
def ivw_from_h(h, min_snp=2):
    if len(h) < min_snp: return None
    flip = h["beta"] < 0
    h = h.copy()
    h.loc[flip, "beta"] = -h.loc[flip, "beta"]
    h.loc[flip, "beta_out_aligned"] = -h.loc[flip, "beta_out_aligned"]
    h["theta"] = h["beta_out_aligned"] / h["beta"]
    w = (h["beta"] ** 2) / (h["se_out"] ** 2)
    theta = float((w * h["theta"]).sum() / w.sum())
    se = float(np.sqrt(1 / w.sum()))
    return dict(beta_MR=theta, se=se, or_=np.exp(theta),
                or_lo=np.exp(theta - 1.96 * se), or_hi=np.exp(theta + 1.96 * se),
                p=z2p(theta / se))

adult_rows = []
# GLP1R: adult has NO own instrument (top p=3.2e-7); fixed infancy union looked up in adult BMI
# (same logic as 7y/8y). Expect weak-instrument explosion -> report numbers but mark uninterpretable.
awin = pd.read_csv(f"{DIR}/AdultBMI_ieub40_GLP1R_window.csv")
minp_a = float(awin[awin["in_100kb"]]["p"].min())
inst_g = pd.concat([pd.read_csv(f) for f in glob.glob(f"{DIR}/instruments_MoBa2022_*_GLP1R_100kb.csv")])
union_g = inst_g.drop_duplicates("rsid")[["rsid"]]
exp_a = awin[awin["rsid"].isin(union_g["rsid"])].copy()
meanF_a = float(((exp_a["beta"] / exp_a["se"]) ** 2).mean())
for model in ["Model1", "Model2", "Model3"]:
    outc = pd.read_csv(f"{DIR}/iPSYCH_{model}_GLP1R_window.csv")
    h = harmonize(exp_a, outc)
    r = ivw_from_h(h)
    adult_rows.append(dict(gene="GLP1R", timepoint="adult", age=18, model=model,
                           n_snp=len(h), mean_F=meanF_a, min_p_100kb=minp_a, **(r or {})))
print(f"[adult GLP1R] fixed union lookup: {len(exp_a)} SNPs found, mean F={meanF_a:.2f} (weak)")
# MC4R: own adult instruments (already harmonized mrinput) -- kept in CSV for archive only
for model in ["Model1", "Model2", "Model3"]:
    h = pd.read_csv(f"{DIR}/mrinput_AdultBMI_ieub40_MC4R_100kb_vs_iPSYCH_{model}.csv")
    meanF_m = float(((h["beta"] / h["se"]) ** 2).mean())
    r = ivw_from_h(h)
    adult_rows.append(dict(gene="MC4R", timepoint="adult", age=18, model=model,
                           n_snp=len(h), mean_F=meanF_m, **(r or {})))
adult_df = pd.DataFrame(adult_rows)
print("\n=== adult (ieu-b-40) ===")
print(adult_df[["gene", "model", "n_snp", "or_", "or_lo", "or_hi", "p", "mean_F"]].to_string(index=False))
res = pd.concat([res, adult_df], ignore_index=True)
# reciprocal: drug-mimicking direction, per -1 SD BMI
res["or_down"] = 1 / res["or_"]
res["or_down_lo"] = 1 / res["or_hi"]
res["or_down_hi"] = 1 / res["or_lo"]
res.to_csv(f"{DIR}/child_mr_developmental.csv", index=False)
res.to_csv(f"{WS}/child_mr_developmental.csv", index=False)
pd.set_option("display.width", 220)
for gene in ["GLP1R", "MC4R"]:
    print(f"\n=== {gene} Model1 ===")
    print(res[(res["gene"] == gene) & (res["model"] == "Model1")]
          [["timepoint", "n_snp", "or_", "or_lo", "or_hi", "p", "mean_F"]].to_string(index=False))

# ---- trajectory plot: GLP1R only, publication layout (2026-09-27 v3) ----
# v3 fixes: adult exploded ORs no longer plotted as clipped triangles (was 3 overlapping
# annotations); legend moved outside axes; adult slot gets one tidy note box.
from matplotlib.lines import Line2D
plt.rcParams["font.family"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
colors = {"Model1": "#2166ac", "Model2": "#4dac26", "Model3": "#d01c8b"}
labels = {"Model1": "M1 主分析", "Model2": "M2 校正精神科诊断", "Model3": "M3 情感障碍亚组"}
DODGE = {"Model1": -0.09, "Model2": 0.0, "Model3": 0.09}
X_META = 18.0
XT = [0, 1, 2, 3, 5, 7, 8, X_META]
XTL = ["0", "1", "2", "3", "5", "7", "8", "成人\n(meta)"]

# adult layer: own instruments from BMI_meta (ieu-b-40+MVP meta, N~1.1M); MVP = sensitivity only
adult_mm = pd.read_csv(f"{DIR}/adult_meta_mvp_mr.csv")
adult_mm = adult_mm[adult_mm["dataset"] == "AdultBMI_meta"]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.6), gridspec_kw={"width_ratios": [1.55, 1]})
g = "GLP1R"
for a in (ax, ax2):
    a.axvspan(0.15, 3.5, color="#2166ac", alpha=0.05, zorder=0)   # valid window (F>=10)
    a.axvline(13, color="lightgrey", lw=1, zorder=0)              # separator before adult

handles = {}
for model in ["Model1", "Model2", "Model3"]:
    sub = res[(res["gene"] == g) & (res["model"] == model) & (res["timepoint"] != "adult")] \
        .dropna(subset=["or_"]).sort_values("age")
    x = sub["age"] + DODGE[model]
    strong = sub["mean_F"] >= 10
    s_s, s_w = sub[strong], sub[~strong]
    handles[labels[model]] = ax.errorbar(
        x[strong], s_s["or_down"],
        yerr=[s_s["or_down"] - s_s["or_down_lo"], s_s["or_down_hi"] - s_s["or_down"]],
        fmt="o-", color=colors[model], ecolor=colors[model], elinewidth=1.2,
        capsize=2.5, ms=4.5, lw=1.2, alpha=0.95)
    if len(s_w):
        ax.errorbar(x[~strong], s_w["or_down"],
                    yerr=[s_w["or_down"] - s_w["or_down_lo"], s_w["or_down_hi"] - s_w["or_down"]],
                    fmt="o", mfc="none", mec=colors[model], ecolor=colors[model],
                    elinewidth=1.0, capsize=2.5, ms=5, lw=0, alpha=0.45)
handles["F<10 弱工具（不可解释）"] = Line2D([0], [0], marker="o", mfc="none", mec="grey",
                                            color="grey", ls="none", ms=5)

# adult points: meta squares, model colours, own instruments (F=44)
sub = adult_mm.set_index("model")
for model in ["Model1", "Model2", "Model3"]:
    if model not in sub.index: continue
    r = sub.loc[model]
    ax.errorbar(X_META + DODGE[model] * 2.2, r["or_down"],
                yerr=[[r["or_down"] - r["or_down_lo"]], [r["or_down_hi"] - r["or_down"]]],
                fmt="s", color=colors[model], ecolor=colors[model], elinewidth=1.3,
                capsize=3, ms=6.5, lw=0, alpha=0.95)
    if r["p"] < 0.05:
        ax.text(X_META + DODGE[model] * 2.2, r["or_down_hi"] * 1.25, "*", ha="center",
                fontsize=13, color=colors[model], fontweight="bold")
handles["成人 BMI_meta（N=110万，自有工具）"] = Line2D([0], [0], marker="s", color="grey", ls="none", ms=6.5)

ax.axhline(1, color="grey", ls="--", lw=1)
ax.set_yscale("log")
ax.set_ylim(0.09, 20)
ax.set_yticks([0.125, 0.25, 0.5, 1, 2, 4, 8, 16])
ax.set_yticklabels(["0.125", "0.25", "0.5", "1", "2", "4", "8", "16"])
ax.set_xticks(XT); ax.set_xticklabels(XTL)
ax.set_xlim(-0.8, 20.5)
ax.set_xlabel("儿童年龄（岁）")
ax.set_ylabel("OR（每降低 1 SD BMI → 自杀未遂）")
ax.set_title("A  GLP1R 固定工具 × 发育时间窗（药效方向）", fontsize=12, loc="left")
ax.text(X_META, 0.30, "成人层 = BMI_meta（N≈110万，自有工具）",
        ha="center", va="center", fontsize=8, color="#555555",
        bbox=dict(boxstyle="round,pad=0.45", fc="#f2f2f2", ec="lightgrey"))
ax.text(1.82, 15, "有效工具窗（F≥10）", fontsize=8.5, color="#2166ac", ha="center", alpha=0.85)

m1 = res[(res["gene"] == g) & (res["model"] == "Model1") & (res["timepoint"] != "adult")].sort_values("age")
ax2.plot(m1["age"], m1["mean_F"], "o-", color="#b2182b", lw=1.5, ms=5)
frow = adult_mm[adult_mm["model"] == "Model1"]
if len(frow):
    fa = float(frow["mean_F"].iloc[0])
    ax2.plot(X_META, fa, "s", color="#b2182b", mec="black", ms=8)
    ax2.annotate(f"成人 meta F={fa:.0f}", xy=(X_META, fa), xytext=(X_META - 4.2, fa * 1.4),
                 fontsize=8.5, color="#b2182b",
                 arrowprops=dict(arrowstyle="-", color="grey", lw=0.8))
ax2.axhline(10, color="grey", ls=":", lw=1)
ax2.text(5.2, 11.3, "F=10 弱工具阈值", fontsize=8.5, color="grey")
ax2.set_yscale("log"); ax2.set_yticks([1, 3, 10, 30, 100]); ax2.set_yticklabels(["1", "3", "10", "30", "100"])
ax2.set_xticks(XT); ax2.set_xticklabels(XTL)
ax2.set_xlim(-0.8, 20.5)
ax2.set_xlabel("儿童年龄（岁）"); ax2.set_ylabel("工具平均 F 统计量")
ax2.set_title("B  工具强度随年龄变化", fontsize=12, loc="left")

fig.suptitle("GLP1R 顺式工具 × 儿童期 BMI 发育时间窗 → iPSYCH 自杀未遂（<32 岁）", fontsize=13.5, y=0.995)
fig.legend(handles.values(), handles.keys(), loc="upper center", ncol=3, fontsize=9,
           frameon=False, bbox_to_anchor=(0.5, 0.945))
fig.tight_layout(rect=[0, 0, 1, 0.855])
for out in [f"{DIR}/child_mr_developmental.png", f"{WS}/child_mr_developmental.png"]:
    fig.savefig(out, dpi=300, bbox_inches="tight")
print("\nsaved child_mr_developmental.csv / .png")
