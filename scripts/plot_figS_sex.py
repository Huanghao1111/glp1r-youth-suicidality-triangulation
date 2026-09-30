#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Supplementary figure: sex-stratified YRBS 2025 associations (v2, clean layout).
Panel A: trying-to-lose-weight -> 3 suicidality outcomes by sex (all BMI + normal-BMI strata).
Panel B: perceived overweight at normal BMI -> outcomes, by sex.
Output: results/FigS_sex_stratified.png/.pdf (300 dpi)."""
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

SEXC = {"Female": "#B35806", "Male": "#2166AC"}
OUTLAB = {"persistent_sadness": "Persistent sadness /\nhopelessness",
          "seriously_considered": "Seriously considered\nsuicide",
          "suicide_attempt": "Suicide attempt"}
OUTORDER = ["persistent_sadness", "seriously_considered", "suicide_attempt"]

df = pd.read_csv("results/yrbs_2025_sex_interaction.csv")

def fmt_p(p):
    return "<.001" if p < 0.001 else f"={p:.2f}".replace("0.", ".")

def get(section, outcome, sex=None, stratum=None):
    q = df[(df.section == section) & (df.outcome == outcome)]
    if sex is not None:
        q = q[q.sex == sex]
    if stratum is not None:
        q = q[q.stratum == stratum]
    return q.iloc[0] if len(q) else None

fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 4.1),
                               gridspec_kw={"width_ratios": [1.5, 1.0], "wspace": 0.55})

# ---------- Panel A: build row list top-to-bottom ----------
rowsA = []
for oc in OUTORDER:
    rowsA.append(("header", OUTLAB[oc], None))
    for scope, slab in [("all_BMI", "all BMI"), ("Normal", "normal BMI")]:
        for sx in ["Female", "Male"]:
            rowsA.append(("row", f"   {sx}, {slab}", get("stratified", oc, sx, scope)))
    p = get("interaction", oc, stratum="all_BMI").p
    rowsA.append(("pint", f"P(interaction){fmt_p(p)}", None))
    rowsA.append(("gap", "", None))

nA = len(rowsA)
for i, (kind, lab, r) in enumerate(rowsA):
    y = nA - i
    if kind == "row" and r is not None:
        c = SEXC[lab.strip().split(",")[0]]
        normal = "normal BMI" in lab
        axA.errorbar(r.OR, y, xerr=[[r.OR - r.lo], [r.hi - r.OR]], fmt="o",
                     ms=4.2 if normal else 3.8,
                     mfc=c if normal else "white", mec=c, ecolor=c,
                     elinewidth=0.9, capsize=1.8, capthick=0.9, zorder=3)
        axA.text(8.4, y, f"{r.OR:.2f} ({r.lo:.2f}\u2013{r.hi:.2f})",
                 va="center", ha="right", fontsize=6.9,
                 color="#111111" if r.p < 0.05 else "#666666")

axA.axvline(1.0, color="#888888", lw=0.7, ls="--", zorder=1)
axA.set_xscale("log"); axA.set_xlim(0.75, 9.2)
axA.set_xticks([1, 2, 4, 8]); axA.set_xticklabels(["1", "2", "4", "8"])
axA.set_yticks([nA - i for i in range(nA)])
axA.set_yticklabels([lab for _, lab, _ in rowsA], fontsize=6.9)
for i, (kind, lab, _) in enumerate(rowsA):
    tl = axA.get_yticklabels()[i]
    if kind == "header":
        tl.set_fontweight("bold"); tl.set_fontsize(7.2)
    elif kind == "pint":
        tl.set_fontstyle("italic"); tl.set_fontsize(6.3); tl.set_color("#555555")
axA.set_ylim(0.3, nA + 0.7)
axA.set_xlabel("Adjusted OR (95% CI), trying to lose weight vs not", fontsize=7.3)
axA.set_title("A  Weight-loss intention \u2192 suicidality, by sex", fontsize=7.8,
              fontweight="bold", loc="left")
axA.spines[["top", "right"]].set_visible(False)
hF = axA.errorbar([], [], yerr=1, fmt="o", ms=4.0, mec=SEXC["Female"], ecolor=SEXC["Female"], mfc=SEXC["Female"], label="Female")
hM = axA.errorbar([], [], yerr=1, fmt="o", ms=4.0, mec=SEXC["Male"], ecolor=SEXC["Male"], mfc=SEXC["Male"], label="Male")
hO = axA.errorbar([], [], yerr=1, fmt="o", ms=3.8, mec="#333333", ecolor="#333333", mfc="white", label="all BMI (open)")
hN = axA.errorbar([], [], yerr=1, fmt="o", ms=4.2, mec="#333333", ecolor="#333333", mfc="#333333", label="normal BMI (filled)")
axA.legend(handles=[hF, hM, hO, hN], loc="lower left", bbox_to_anchor=(0.005, 0.03),
           frameon=False, fontsize=6.2, handletextpad=0.15, borderaxespad=0.0,
           labelspacing=0.25, ncol=2, columnspacing=0.6)

# ---------- Panel B ----------
rowsB = []
for oc in OUTORDER:
    rowsB.append(("header", OUTLAB[oc], None))
    for sx in ["Female", "Male"]:
        rowsB.append(("row", f"   {sx}", get("discordance_normalBMI", oc, sx, "Normal")))
    p = get("interaction_discordance", oc, stratum="Normal_only").p
    rowsB.append(("pint", f"P(interaction){fmt_p(p)}", None))
    rowsB.append(("gap", "", None))

nB = len(rowsB)
for i, (kind, lab, r) in enumerate(rowsB):
    y = nB - i
    if kind == "row" and r is not None:
        c = SEXC[lab.strip()]
        axB.errorbar(r.OR, y, xerr=[[r.OR - r.lo], [r.hi - r.OR]], fmt="s",
                     ms=4.0, color=c, ecolor=c, elinewidth=0.9,
                     capsize=1.8, capthick=0.9, zorder=3)
        axB.text(min(r.hi * 1.22, 26), y, f"{r.OR:.2f} ({r.lo:.2f}\u2013{r.hi:.2f})",
                 va="center", ha="left", fontsize=6.9,
                 color="#111111" if r.p < 0.05 else "#666666")

axB.axvline(1.0, color="#888888", lw=0.7, ls="--", zorder=1)
axB.set_xscale("log"); axB.set_xlim(0.75, 32)
axB.set_xticks([1, 2, 4, 8]); axB.set_xticklabels(["1", "2", "4", "8"])
axB.set_yticks([nB - i for i in range(nB)])
axB.set_yticklabels([lab for _, lab, _ in rowsB], fontsize=6.9)
for i, (kind, lab, _) in enumerate(rowsB):
    tl = axB.get_yticklabels()[i]
    if kind == "header":
        tl.set_fontweight("bold"); tl.set_fontsize(7.2)
    elif kind == "pint":
        tl.set_fontstyle("italic"); tl.set_fontsize(6.3); tl.set_color("#555555")
axB.set_ylim(0.3, nB + 0.7)
axB.set_xlabel("Adjusted OR (95% CI)", fontsize=7.3)
axB.set_title("B  Perceived overweight at normal BMI, by sex", fontsize=7.8,
              fontweight="bold", loc="left")
axB.spines[["top", "right"]].set_visible(False)

fig.savefig("results/FigS_sex_stratified.png", dpi=300, bbox_inches="tight")
fig.savefig("results/FigS_sex_stratified.pdf", bbox_inches="tight")
print("saved results/FigS_sex_stratified.png/.pdf")
