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

# ---------------- Figure 1: FAERS ----------------

by = pd.read_csv("results/faers_full_ped_glp1_byyear.csv")
rp = pd.read_csv("results/faers_full_ror_pediatric.csv")
ra = pd.read_csv("results/faers_full_ror_adult.csv")
rp = rp[rp["drug"] == "GLP1-class(all)"].set_index("event_group")
ra = ra[ra["drug"] == "GLP1-class(all)"].set_index("event_group")
groups = ["Suicidality", "Depression", "Anxiety", "Eating disorder", "Psychosis"]
glab = {"Suicidality": "Suicidality", "Depression": "Depression", "Anxiety": "Anxiety",
        "Eating disorder": "Eating disorders", "Psychosis": "Psychotic symptoms"}

drugs = [("sema", "Semaglutide", "#C2185B"), ("lira", "Liraglutide", "#E65100"),
         ("tirz", "Tirzepatide", "#7B1FA2"), ("setm", "Setmelanotide", "#00897B"),
         ("dula", "Dulaglutide", "#3D5AFE"), ("exen", "Exenatide", "#558B2F"),
         ("lixi", "Lixisenatide", "#795548")]

fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.9, 3.3), gridspec_kw={"width_ratios": [1.05, 1.35]})

# Panel A: stacked annual bars
years = by["yr"].astype(int).tolist()
bottoms = np.zeros(len(by))
for col, name, color in drugs:
    vals = by[col].fillna(0).values
    a1.bar(years, vals, bottom=bottoms, color=color, width=0.82, label=name,
           edgecolor="white", linewidth=0.3)
    bottoms += vals
for x, tot in zip(years, bottoms):
    if tot > 0:
        a1.text(x, tot + 2.5, str(int(tot)), ha="center", va="bottom", fontsize=6.3, color="#333333")
a1.set_xticks(years)
a1.set_xticklabels([str(yy) for yy in years[:-1]] + ["2026*"], rotation=45, ha="right", fontsize=6.8)
a1.set_ylabel("Pediatric GLP-1 RA reports, n", fontsize=7.6)
a1.set_ylim(0, max(bottoms) * 1.14)
for s in ["top", "right"]:
    a1.spines[s].set_visible(False)
a1.legend(fontsize=5.9, frameon=False, loc="upper left", handlelength=1.0,
          handletextpad=0.4, labelspacing=0.35, borderaxespad=0.0)
a1.set_title("A", loc="left", fontweight="bold", fontsize=9)

# Panel B: ROR forest ped vs adult
y = 0.0
for g in reversed(groups):
    y += 1.0
    for dfc, mk, col, dy in [(ra, "s", "#1565C0", -0.16), (rp, "o", "#C2185B", 0.16)]:
        r = dfc.loc[g]
        a2.errorbar(r["ROR"], y + dy,
                    xerr=[[r["ROR"] - r["ROR_lo"]], [r["ROR_hi"] - r["ROR"]]],
                    fmt=mk, ms=4.0, color=col, ecolor=col, elinewidth=1.0,
                    capsize=2.0, capthick=1.0, zorder=3)
        star = "*" if r["ROR_lo"] > 1 or r["ROR_hi"] < 1 else ""
        a2.text(23, y + dy, f"a={int(r['a'])}   {r['ROR']:.2f} ({r['ROR_lo']:.2f}\u2013{r['ROR_hi']:.2f}){star}",
                va="center", ha="left", fontsize=6.4, color=col)
a2.axvline(1.0, ls=(0, (4, 3)), lw=0.8, color="#555555", zorder=1)
a2.set_xscale("log")
a2.set_xlim(0.045, 22)
a2.set_ylim(0.4, y + 0.7)
a2.set_yticks([i + 0 for i in range(1, len(groups) + 1)])
a2.set_yticklabels([glab[g] for g in reversed(groups)], fontsize=7.4)
a2.set_xticks([0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20])
a2.set_xticklabels(["0.05", "0.1", "0.2", "0.5", "1", "2", "5", "10", "20"], fontsize=6.8)
for s in ["top", "right"]:
    a2.spines[s].set_visible(False)
a2.set_xlabel("ROR (95% CI), log scale;  * 95% CI excludes 1", fontsize=7.6)
a2.set_title("B", loc="left", fontweight="bold", fontsize=9)
handles = [plt.Line2D([], [], marker="o", ls="", color="#C2185B", ms=4.5,
                      label="Pediatric (<18 y; N=616,860)"),
           plt.Line2D([], [], marker="s", ls="", color="#1565C0", ms=4.2,
                      label="Adult (\u226518 y; N=9,457,685)")]
a2.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, -0.42), fontsize=6.2,
          frameon=False, ncol=2, columnspacing=1.0, handletextpad=0.3)
fig.tight_layout(w_pad=2.2, rect=[0, 0.02, 1, 1])
import os
os.makedirs("results/pub_final", exist_ok=True)
fig.savefig("results/pub_final/Figure_1.png", dpi=600)
fig.savefig("results/pub_final/Figure_1.pdf")
plt.close(fig)
print("Fig3 saved")
