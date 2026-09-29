#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""JAACAP graphical abstract: triangulation of GLP-1 RA youth suicidality evidence.
Output: results/pub_final/Graphical_Abstract.png (2656x1062 @300dpi) + .pdf"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
import numpy as np

for f in ["C:/Windows/Fonts/arial.ttf", "C:/Windows/Fonts/arialbd.ttf", "C:/Windows/Fonts/ariali.ttf"]:
    try:
        fm.fontManager.addfont(f)
    except Exception:
        pass
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
                     "text.color": "#1a1a1a"})

C_FAERS = "#D81B8C"   # magenta
C_YRBS  = "#2166AC"   # blue
C_MR    = "#66A61E"   # green
C_DARK  = "#1a1a1a"
C_GREY  = "#555555"

fig = plt.figure(figsize=(8.853, 3.54), dpi=300)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 40); ax.axis("off")

def box(x, y, w, h, fc, ec, lw=1.2, r=1.2):
    b = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                       fc=fc, ec=ec, lw=lw, mutation_aspect=1.0)
    ax.add_patch(b); return b

def arrow(x1, y1, x2, y2, color, lw=2.4, style="-|>", ms=16):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                 mutation_scale=ms, lw=lw, color=color, zorder=5))

# ---------- Title strip ----------
ax.text(50, 37.6, "GLP-1 Receptor Agonists and Suicide-Attempt Risk in Youth: Triangulation of Three Evidence Layers",
        ha="center", va="center", fontsize=10.5, fontweight="bold", color=C_DARK)

# ---------- Left column: three evidence layers ----------
layers = [
    ("FAERS", "Pediatric reports", "2012–2026, 55 quarters",
     "Suicidality signal in children: ROR 1.74\n(1.11–2.72; IC025<0); adults 0.40 (0.38–0.42)\nEmerged post-2023, exceeding the ≤1.5-fold\nrise in stable-warning comparators", C_FAERS, 26.2, 6.6, 2.75),
    ("YRBS 2025", "US adolescents", "survey-weighted",
     "Trying to lose weight:\nsuicide attempt OR 2.28 (1.75–2.97)\nacross BMI categories", C_YRBS, 16.05, 6.9, 3.1),
    ("GLP1R cis-MR", "Drug-target genetics", "MoBa × iPSYCH × GIANT/MVP",
     "BMI-lowering GLP1R action in infancy:\nsuicide attempt OR ≈1.8–2.0 per 1-SD lower BMI\nAdult estimates attenuated, imprecise", C_MR, 5.9, 6.9, 3.1),
]
LW = 26
for tag, sub, src, txt, c, y, fs, dy in layers:
    box(2, y, LW, 9.4, "#FFFFFF", c, lw=1.6)
    box(2, y+7.1, LW, 2.3, c, c, lw=0)
    ax.text(2.9, y+8.25, tag, fontsize=9.5, fontweight="bold", color="white", va="center")
    ax.text(14.6, y+8.25, sub, fontsize=6.6, color="white", va="center", style="italic")
    ax.text(2.9, y+5.9, src, fontsize=6.4, color=C_GREY, va="center")
    ax.text(2.9, y+dy, txt, fontsize=fs, color=C_DARK, va="center", linespacing=1.25)
    arrow(2+LW+0.5, y+4.7, 38.5, y+4.7, c, lw=2.2)

# ---------- Middle: convergence ----------
box(39.5, 12.0, 18.5, 16.5, "#F5F5F7", "#888888", lw=1.2)
ax.text(48.75, 26.0, "TRIANGULATION", ha="center", fontsize=8.6, fontweight="bold", color=C_DARK)
ax.text(48.75, 23.6, "largely non-overlapping biases", ha="center", fontsize=6.4, color=C_GREY, style="italic")
# simple triad icon
tri = [(45.9, 16.6, C_FAERS), (51.6, 16.6, C_YRBS), (48.75, 19.9, C_MR)]
for cx, cy, c in tri:
    ax.add_patch(Circle((cx, cy), 1.45, fc=c, ec="none", zorder=4))
ax.plot([45.9, 51.6], [16.6, 16.6], color="#999999", lw=1.0, zorder=3)
ax.plot([45.9, 48.75], [16.6, 19.9], color="#999999", lw=1.0, zorder=3)
ax.plot([51.6, 48.75], [16.6, 19.9], color="#999999", lw=1.0, zorder=3)
ax.add_patch(Circle((48.75, 17.9), 0.62, fc=C_DARK, ec="none", zorder=5))
ax.text(48.75, 13.4, "directionally convergent signals in youth", ha="center", fontsize=6.4, color=C_DARK)
arrow(58.4, 20.2, 61.8, 20.2, C_DARK, lw=3.0, ms=20)

# ---------- Right: core finding ----------
box(62.5, 12.0, 35.5, 16.5, "#FFFFFF", C_DARK, lw=1.6)
ax.text(64.0, 26.2, "Core finding", fontsize=9.5, fontweight="bold", color=C_DARK)
ax.text(64.0, 21.6,
        "Genetically modeled GLP1R agonism that lowers BMI\nraises suicide-attempt risk in early development\n(infant waves; adult estimates attenuated, imprecise)",
        fontsize=7.4, va="center", linespacing=1.5)
# mini developmental curve
mx0, my0, mw, mh = 65.0, 13.6, 15.5, 4.2
ages = np.array([0, 0.25, 0.5, 0.67, 1.0, 1.5])
ors  = np.array([1.9, 1.81, 1.85, 1.97, 1.93, 1.97])
xsc = mx0 + (ages / 1.5) * mw * 0.62
ysc = my0 + (np.log2(ors) / 1.2) * mh
ax.plot(np.append(xsc, mx0+mw*0.78), np.append(ysc, my0+0.25*mh), "-o", ms=3.4, lw=1.6,
        color=C_MR, mfc="white", mew=1.0, zorder=6)
ax.plot([mx0, mx0+mw], [my0, my0], ls=(0, (4, 3)), lw=0.8, color="#888888")
ax.annotate("adult: attenuated, n.s.", xy=(mx0+mw*0.78, my0+0.25*mh),
            xytext=(mx0+mw*0.30, my0-1.35), fontsize=5.8, color=C_GREY,
            arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.7))
ax.text(mx0-1.2, my0+mh*0.55, "OR", fontsize=6.0, rotation=90, color=C_GREY, va="center")
ax.text(mx0+mw*0.28, my0+mh+0.35, "infancy → 1.5 y", fontsize=5.8, color=C_GREY, ha="center")
ax.text(83.2, 16.2, "Early life, per 1-SD\nlower BMI: suicide-attempt\nrisk ≈1.8–2.0× higher", fontsize=6.8,
        fontweight="bold", color=C_MR, va="center", linespacing=1.4)

# ---------- Bottom implication strip ----------
box(2, 0.8, 96, 4.6, "#1a1a1a", "#1a1a1a", lw=0)
ax.text(50, 3.35, "Clinical implication:  suicide-risk monitoring should accompany GLP-1 RA treatment in adolescents,",
        ha="center", fontsize=8.2, color="white", fontweight="bold")
ax.text(50, 1.55, "and pediatric-specific safety studies are needed; current regulatory evidence is almost entirely adult.",
        ha="center", fontsize=7.4, color="#CCCCCC")

import os
os.makedirs("results/pub_final", exist_ok=True)
fig.savefig("results/pub_final/Graphical_Abstract.png", dpi=300)
fig.savefig("results/pub_final/Graphical_Abstract.pdf")
print("GA saved")
