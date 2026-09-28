#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Publication-grade English figures for the GLP1R youth suicidality manuscript.
Fig 2: wave-specific cis-MR forest (MoBa waves x iPSYCH Models 1-3), OR per +1 SD BMI.
Fig 3: FAERS pediatric pharmacovigilance (A: annual reports by drug; B: ROR ped vs adult).
Output: results/mr_cis/Fig2_forest_pub.png/pdf, results/Fig3_faers_pub.png/pdf (300 dpi)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import pandas as pd

for f in ["C:/Windows/Fonts/arial.ttf"]:
    try:
        fm.fontManager.addfont(f)
    except Exception:
        pass
plt.rcParams.update({
    "font.family": "Arial", "font.size": 7.5, "axes.linewidth": 0.7,
    "axes.edgecolor": "#222222", "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.color": "#222222", "ytick.color": "#222222", "text.color": "#111111",
    "axes.labelcolor": "#111111", "pdf.fonttype": 42, "ps.fonttype": 42,
})

MC = {"Model1": "#2166AC", "Model2": "#66A61E", "Model3": "#D81B8C"}
MLAB = {"Model1": "Model 1 (primary)", "Model2": "Model 2 (adj. psychiatric dx)",
        "Model3": "Model 3 (affective subgroup)"}

# ---------------- Figure S1: cis-MR forest (inverted) ----------------
df = pd.read_csv("results/mr_cis/child_mr_results.csv")
df = df[(df["gene"] == "GLP1R") & (df["window"] == "100kb") & df["dataset"].str.startswith("MoBa2022")]
order = ["3months", "6months", "8months", "1year", "1.5years"]
wlab = {"3months": "3 months", "6months": "6 months", "8months": "8 months",
        "1year": "1 year", "1.5years": "1.5 years"}
df["wave"] = df["dataset"].str.replace("MoBa2022_", "", regex=False)
df = df[df["wave"].isin(order)]

rows = []
for w in order:  # bottom-to-top: 1.5y first (top), so reverse later
    for m in ["Model1", "Model2", "Model3"]:
        r = df[(df["wave"] == w) & (df["model"] == m)]
        if len(r):
            rows.append((w, m, r.iloc[0]))

fig, ax = plt.subplots(figsize=(7.6, 4.6))
ypos, labels, ann = [], [], []
y = 0.0
for w in reversed(order):
    block = [r for r in rows if r[0] == w]
    if not block:
        continue
    for (ww, m, r) in block:
        y += 1.0
        orr, lo, hi = 1.0/r["or_"], 1.0/r["or_hi"], 1.0/r["or_lo"]
        ax.errorbar(orr, y, xerr=[[orr - lo], [hi - orr]],
                    fmt="o", ms=4.2, color=MC[m], ecolor=MC[m], elinewidth=1.0,
                    capsize=2.2, capthick=1.0, zorder=3)
        sig = r["p"] < 0.05
        txt = f"{orr:.2f} ({lo:.2f}\u2013{hi:.2f})"
        ax.text(6.85, y, txt, clip_on=False, va="center", ha="right", fontsize=6.8,
                fontweight="bold" if sig else "normal",
                color="#111111" if sig else "#666666")
        ptxt = f"{r['p']:.3f}" if r["p"] >= 0.001 else "<0.001"
        ax.text(8.2, y, ptxt, clip_on=False,
                va="center", ha="right", fontsize=6.8,
                fontweight="bold" if sig else "normal",
                color="#111111" if sig else "#666666")
        labels.append((y, f"{m.replace('Model', 'M')} (n={int(r['n_snp'])})"))
    ax.text(-0.012, y + 0.62, wlab[w], transform=ax.get_yaxis_transform(),
            rotation=0, va="center", ha="right", fontsize=7.6, fontweight="bold")
    y += 0.9  # gap between waves

for yy, lab in labels:
    ax.text(-0.012, yy, lab, transform=ax.get_yaxis_transform(), va="center", ha="right", fontsize=7.2)

ax.axvline(1.0, ls=(0, (4, 3)), lw=0.8, color="#555555", zorder=1)
ax.set_xscale("log")
ax.set_xlim(0.85, 5.0)
ax.set_ylim(0.2, y + 0.6)
ax.set_xticks([1.0, 1.5, 2.0, 3.0, 5.0])
ax.set_xticklabels(["1.0", "1.5", "2.0", "3.0", "5.0"])
from matplotlib.ticker import NullFormatter
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_yticks([])
for s in ["left", "top", "right"]:
    ax.spines[s].set_visible(False)
ax.set_xlabel("Odds ratio for suicide attempt per 1-SD lower childhood BMI z-score (log scale)", fontsize=7.8)
ax.text(6.62, y + 0.35, "OR (95% CI)", ha="right", fontsize=7.4, fontweight="bold", clip_on=False)
ax.text(8.28, y + 0.35, "P value", ha="right", fontsize=7.4, fontweight="bold", clip_on=False)

handles = [plt.Line2D([], [], marker="o", ls="", color=MC[m], ms=4.5, label=MLAB[m]) for m in MC]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.45, 0.0), fontsize=6.8,
           frameon=False, ncol=3, handletextpad=0.2, columnspacing=1.2)
fig.tight_layout(rect=[0, 0.055, 0.78, 1])
fig.savefig("results/pub_final/Figure_S1.png", dpi=600)
fig.savefig("results/pub_final/Figure_S1.pdf")
plt.close(fig)
print("Fig2 saved")

