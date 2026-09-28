#!/usr/bin/env python3
"""Forest plot: FAERS pediatric vs adult RORs."""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.sans-serif"] = ["Noto Sans CJK SC", "SimHei", "Microsoft YaHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

def load(path):
    return list(csv.DictReader(open(path)))

ped = load("results/faers_ror_pediatric.csv")
adu = load("results/faers_ror_adult.csv")

groups = ["Suicidality", "Depression", "Anxiety", "Eating disorder", "Psychosis"]
zh = {"Suicidality": "自杀相关", "Depression": "抑郁", "Anxiety": "焦虑",
      "Eating disorder": "进食障碍", "Psychosis": "精神病性症状"}

def get(rows, drug, grp):
    for r in rows:
        if r["drug"] == drug and r["event_group"] == grp:
            return float(r["ror"]), float(r["lo"]), float(r["hi"]), int(r["a"])
    return None

fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), gridspec_kw={"width_ratios": [1, 1.35]})

# Panel A: GLP-1 class, pediatric vs adult
ax = axes[0]
y = np.arange(len(groups))[::-1] * 1.0
for i, (rows, lab, col, off) in enumerate([(ped, "儿科 (<18 岁)", "#C0392B", 0.14),
                                            (adu, "成人 (≥18 岁)", "#2E5E8C", -0.14)]):
    for yi, grp in zip(y, groups):
        g = get(rows, "GLP1-class(all)", grp)
        if not g: continue
        v, lo, hi, a = g
        ax.errorbar(v, yi + off, xerr=[[v - lo], [hi - v]], fmt="o", color=col,
                    ecolor=col, elinewidth=1.6, capsize=3, ms=6,
                    label=lab if yi == y[0] else None)
        ax.annotate(f"n={a}", (hi * 1.15, yi + off), va="center", fontsize=8, color=col)
ax.axvline(1, ls="--", c="grey", lw=1)
ax.set_xscale("log")
ax.set_xlim(0.08, 60)
ax.set_xticks([0.1, 1, 10])
ax.set_xticklabels(["0.1", "1", "10"])
ax.set_yticks(y)
ax.set_yticklabels([zh[g] for g in groups])
ax.set_xlabel("ROR (log 尺度, 95% CI)")
ax.set_title("A. GLP-1 类药物整体：儿科 vs 成人", fontsize=11)
ax.legend(loc="lower right", fontsize=9)

# Panel B: pediatric per-drug, cells with a >= 1
ax = axes[1]
sel = []
for drug in ["Semaglutide", "Liraglutide", "Tirzepatide", "Setmelanotide"]:
    for grp in groups:
        g = get(ped, drug, grp)
        if g and g[3] >= 1:
            sel.append((f"{drug} × {zh[grp]}", *g))
sel.sort(key=lambda x: -x[1])
y2 = np.arange(len(sel))[::-1]
for yi, (lab, v, lo, hi, a) in zip(y2, sel):
    sig = lo > 1
    ax.errorbar(v, yi, xerr=[[v - lo], [hi - v]], fmt="s",
                color="#C0392B" if sig else "#888", ecolor="#C0392B" if sig else "#888",
                elinewidth=1.6, capsize=3, ms=6)
    ax.annotate(f"n={a}", (hi * 1.15, yi), va="center", fontsize=8, color="#555")
ax.axvline(1, ls="--", c="grey", lw=1)
ax.set_xscale("log")
ax.set_xlim(0.05, 400)
ax.set_xticks([0.1, 1, 10, 100])
ax.set_xticklabels(["0.1", "1", "10", "100"])
ax.set_yticks(y2)
ax.set_yticklabels([s[0] for s in sel], fontsize=9)
ax.set_xlabel("ROR (log 尺度, 95% CI)")
ax.set_title("B. 儿科分层：单药信号（仅显示报告数 ≥1 的格子）", fontsize=11)

fig.suptitle("FAERS/openFDA 不相称分析（预分析）：减重/降糖药物 × 精神科不良事件", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("results/faers_forest.png", dpi=200)
print("saved -> results/faers_forest.png")
