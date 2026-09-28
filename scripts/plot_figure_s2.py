#!/usr/bin/env python3
"""Quarter-resolved pediatric GLP-1 RA suicidality disproportionality figure
with regulatory milestones. Output: results/pub_final/Figure_S2.png/.pdf
NOTE: this figure appears as Figure S1 in the published supplement
(renumbered at final QC; script/output filenames predate the renumbering)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import pandas as pd

for f in ["C:/Windows/Fonts/arial.ttf", "C:/Windows/Fonts/arialbd.ttf"]:
    try:
        fm.fontManager.addfont(f)
    except Exception:
        pass
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42})

df = pd.read_csv("results/faers_ped_suic_quarterly.csv")
df["t"] = df["quarter"].str.extract(r"(\d{4})Q(\d)").astype(int).apply(lambda r: r[0] + (r[1]-1)/4, axis=1)

C_BAR = "#D81B8C"
C_LINE = "#555555"
C_ROR = "#2166AC"

fig, axes = plt.subplots(2, 1, figsize=(9.2, 5.4), dpi=300,
                         gridspec_kw={"height_ratios": [1, 1.15], "hspace": 0.32},
                         sharex=True)

milestones = [(2022.75, "Semaglutide adolescent\nobesity approval\n(Dec 2022)", 5.35),
              (2023.25, "EMA PRAC review\nannounced\n(Apr 2023)", 4.05),
              (2024.0, "FDA update\n(Jan 2024)", 5.35),
              (2024.25, "EMA PRAC\nconcluded\n(Apr 2024)", 4.05)]

ax = axes[0]
ax.bar(df["t"], df["a"], width=0.22, color=C_BAR, alpha=0.85, label="Pediatric suicidality cases (a)")
ax2 = ax.twinx()
ax2.plot(df["t"], df["n_glp"], color=C_LINE, lw=1.4, label="Pediatric GLP-1 RA reports (n)")
ax2.set_ylabel("Pediatric GLP-1 RA reports per quarter", fontsize=8, color=C_LINE)
ax2.tick_params(axis="y", labelsize=7.5, colors=C_LINE)
ax.set_ylabel("Suicidality cases per quarter", fontsize=8, color=C_BAR)
ax.tick_params(axis="y", labelsize=7.5, colors=C_BAR)
ax.set_ylim(0, 6.5)
for x, lab, ty in milestones:
    for a_ in axes:
        a_.axvline(x, color="#999999", lw=0.7, ls=(0, (3, 3)), zorder=0)
    ax.annotate(lab, xy=(x, 6.3), xytext=(x+0.12, ty), fontsize=6.0, color="#666666")
ax.text(2021.0, 5.6, "A", fontsize=11, fontweight="bold")
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1+h2, l1+l2, fontsize=6.8, loc="upper left", frameon=False,
          bbox_to_anchor=(0.0, 0.82))

ax = axes[1]
ax.axhline(1, color="#333333", lw=0.8)
shown = df[df["t"] >= 2021.0]
for _, r in shown.iterrows():
    fmt = dict(color=C_ROR, lw=1.0) if not r.haldane else dict(color="#999999", lw=0.8)
    ax.plot([r["t"], r["t"]], [r["lo"], min(r["hi"], 400)], **fmt, zorder=2)
    ax.plot(r["t"], r["ROR"], "o", ms=3.6,
            color=C_ROR if not r.haldane else "#999999",
            mfc="white" if r.haldane else C_ROR, mec=C_ROR if not r.haldane else "#999999",
            zorder=3)
ax.set_yscale("log")
ax.set_ylim(0.02, 400)
ax.set_ylabel("Quarterly ROR for suicidality\n(pediatric, class level; log scale)", fontsize=8)
ax.tick_params(axis="y", labelsize=7.5)
ax.tick_params(axis="x", labelsize=7.5)
ax.set_xlim(2020.8, 2026.75)
ax.set_xlabel("FAERS report quarter", fontsize=8.5)
ax.text(2021.0, 210, "B", fontsize=11, fontweight="bold")
ax.text(2021.0, 0.032, "Open grey symbols: Haldane-corrected quarters with zero cases (point estimates not interpretable).\n"
        "Quarters before 2021 had no pediatric GLP-1 RA suicidality cases.",
        fontsize=6.2, color="#666666")

fig.savefig("results/pub_final/Figure_S2.png", dpi=300, bbox_inches="tight")
fig.savefig("results/pub_final/Figure_S2.pdf", bbox_inches="tight")
print("Figure saved")
