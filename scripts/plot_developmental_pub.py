#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Publication-grade developmental figure (JAACAP/JCPP style).
Reads child_mr_developmental.csv + adult_meta_mvp_mr.csv (analysis unchanged).
All-English, Arial, 182 mm double-column width, 600 dpi PNG + vector PDF.
Outputs: E:/CM/GLP1/mr_input/Fig_developmental.{png,pdf} + results/mr_cis/."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

DIR = "E:/CM/GLP1/mr_input"
WS = "results/mr_cis"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.5,
    "axes.linewidth": 0.8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "legend.fontsize": 7,
    "pdf.fonttype": 42,   # editable text in Illustrator
    "ps.fonttype": 42,
})

res = pd.read_csv(f"{DIR}/child_mr_developmental.csv")
adult = pd.read_csv(f"{DIR}/adult_meta_mvp_mr.csv")
adult = adult[adult["dataset"] == "AdultBMI_meta"]

COL = {"Model1": "#2166ac", "Model2": "#4dac26", "Model3": "#d01c8b"}
LAB = {"Model1": "Model 1 (primary)", "Model2": "Model 2 (adj. psych. dx)",
       "Model3": "Model 3 (affective subgr.)"}
DODGE = {"Model1": -0.09, "Model2": 0.0, "Model3": 0.09}
X_ADULT = 18.0
XT = [0, 1, 2, 3, 5, 7, 8, X_ADULT]
XTL = ["0", "1", "2", "3", "5", "7", "8", "Adult"]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.16, 3.45),
                              gridspec_kw={"width_ratios": [1.6, 1], "wspace": 0.28})

for a in (ax, ax2):
    a.axvspan(0.15, 3.5, color="#2166ac", alpha=0.05, zorder=0)
    a.axvline(13, color="0.85", lw=0.7, zorder=0)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)

# ---------------- Panel A ----------------
g = "GLP1R"
for model in ["Model1", "Model2", "Model3"]:
    sub = res[(res["gene"] == g) & (res["model"] == model) & (res["timepoint"] != "adult")] \
        .dropna(subset=["or_"]).sort_values("age")
    x = sub["age"] + DODGE[model]
    strong = sub["mean_F"] >= 10
    s_s, s_w = sub[strong], sub[~strong]
    ax.errorbar(x[strong], s_s["or_down"],
                yerr=[s_s["or_down"] - s_s["or_down_lo"], s_s["or_down_hi"] - s_s["or_down"]],
                fmt="o-", color=COL[model], ecolor=COL[model], elinewidth=0.9,
                capsize=1.8, capthick=0.9, ms=3.2, lw=0.9, mew=0.9, alpha=0.95, zorder=3)
    if len(s_w):
        ax.errorbar(x[~strong], s_w["or_down"],
                    yerr=[s_w["or_down"] - s_w["or_down_lo"], s_w["or_down_hi"] - s_w["or_down"]],
                    fmt="o", mfc="none", mec=COL[model], ecolor=COL[model], elinewidth=0.7,
                    capsize=1.5, capthick=0.7, ms=3.4, lw=0, mew=0.8, alpha=0.45, zorder=2)

# adult meta squares
subA = adult.set_index("model")
for model in ["Model1", "Model2", "Model3"]:
    if model not in subA.index:
        continue
    r = subA.loc[model]
    ax.errorbar(X_ADULT + DODGE[model] * 2.2, r["or_down"],
                yerr=[[r["or_down"] - r["or_down_lo"]], [r["or_down_hi"] - r["or_down"]]],
                fmt="s", color=COL[model], ecolor=COL[model], elinewidth=0.9,
                capsize=1.8, capthick=0.9, ms=4.2, lw=0, mew=0.9, alpha=0.95, zorder=3)

ax.axhline(1, color="0.55", ls="--", lw=0.8, zorder=1)
ax.set_yscale("log")
ax.set_ylim(0.09, 20)
ax.set_yticks([0.125, 0.25, 0.5, 1, 2, 4, 8, 16])
ax.set_yticklabels(["0.125", "0.25", "0.5", "1", "2", "4", "8", "16"])
ax.set_xticks(XT); ax.set_xticklabels(XTL)
ax.set_xlim(-0.8, 20.3)
ax.set_xlabel("Age at BMI measurement (years)")
ax.set_ylabel("OR for suicide attempt\nper 1-SD lower BMI z-score (log scale)")
ax.set_title("A", fontsize=10, fontweight="bold", loc="left", pad=4)
ax.text(0.985, 0.965, "Adult layer: BMI meta-analysis (N≈1.11M)",
        transform=ax.transAxes, ha="right", va="top", fontsize=6.2, color="0.35",
        bbox=dict(boxstyle="round,pad=0.4", fc="0.97", ec="0.85", lw=0.6))
ax.text(1.82, 15.2, "Valid-instrument\nwindow (F≥10)", fontsize=6.2, color="#2166ac",
        ha="center", va="center", alpha=0.9, linespacing=1.3)

# ---------------- Panel B ----------------
m1 = res[(res["gene"] == g) & (res["model"] == "Model1") & (res["timepoint"] != "adult")].sort_values("age")
ax2.plot(m1["age"], m1["mean_F"], "o-", color="#b2182b", lw=1.0, ms=3.4, mew=0.9, zorder=3)
frow = adult[adult["model"] == "Model1"]
if len(frow):
    fa = float(frow["mean_F"].iloc[0])
    ax2.plot(X_ADULT, fa, "s", color="#b2182b", mec="black", mew=0.6, ms=5.5, zorder=3)
    ax2.annotate(f"Adult meta: F={fa:.0f}", xy=(X_ADULT, fa), xytext=(X_ADULT - 5.6, fa * 1.35),
                 fontsize=6.5, color="#b2182b",
                 arrowprops=dict(arrowstyle="-", color="0.5", lw=0.7))
ax2.axhline(10, color="0.55", ls=":", lw=0.8, zorder=1)
ax2.text(4.6, 11.2, "Weak-instrument threshold (F=10)", fontsize=6.2, color="0.4")
ax2.set_yscale("log")
ax2.set_yticks([1, 3, 10, 30, 100]); ax2.set_yticklabels(["1", "3", "10", "30", "100"])
ax2.set_xticks(XT); ax2.set_xticklabels(XTL)
ax2.set_xlim(-0.8, 20.3)
ax2.set_ylim(0.15, 200)
ax2.set_xlabel("Age at BMI measurement (years)")
ax2.set_ylabel("Mean instrument F statistic (log scale)")
ax2.set_title("B", fontsize=10, fontweight="bold", loc="left", pad=4)

# ---------------- shared legend ----------------
# ncol fills column-wise: interleave so rows read [M1 M2 M3] / [weak, adult]
hmap = {
    "M1": Line2D([0], [0], marker="o", color=COL["Model1"], lw=0.9, ms=3.4, mew=0.9,
                 label=LAB["Model1"]),
    "M2": Line2D([0], [0], marker="o", color=COL["Model2"], lw=0.9, ms=3.4, mew=0.9,
                 label=LAB["Model2"]),
    "M3": Line2D([0], [0], marker="o", color=COL["Model3"], lw=0.9, ms=3.4, mew=0.9,
                 label=LAB["Model3"]),
    "weak": Line2D([0], [0], marker="o", mfc="none", mec="0.45", color="0.45", ls="none",
                   ms=3.4, mew=0.8, label="Weak instruments (F<10)"),
    "adult": Line2D([0], [0], marker="s", color="0.45", ls="none", ms=4.2, mew=0.9,
                    label="Adult BMI meta (own instr.)"),
}
# ---------------- legend: inside panel A empty lower band ----------------
handles = [hmap["M1"], hmap["M2"], hmap["M3"], hmap["weak"], hmap["adult"]]
leg = ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.245, 0.015),
                ncol=2, frameon=False, fontsize=6.2, handlelength=1.5,
                columnspacing=0.9, labelspacing=0.45, borderaxespad=0)
fig.tight_layout()
# panel letters above axes (top edge is free; legend sits below)
for a, t in ((ax, "A"), (ax2, "B")):
    a.set_title(t, fontsize=10, fontweight="bold", loc="left", pad=4)

for out in [f"{DIR}/Fig_developmental", f"{WS}/Fig_developmental"]:
    fig.savefig(out + ".png", dpi=600, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
print("saved Fig_developmental.png (600dpi) / .pdf")
