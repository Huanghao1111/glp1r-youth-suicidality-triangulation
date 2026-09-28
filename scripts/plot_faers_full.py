# -*- coding: utf-8 -*-
# FAERS 全量版图: A=儿科GLP-1报告年度趋势(分药堆叠) B=儿科vs成人类级ROR森林图
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.sans-serif"] = ["Noto Sans CJK SC", "SimHei", "Microsoft YaHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

DRUG_CN = {"sema": "Semaglutide", "lira": "Liraglutide", "tirz": "Tirzepatide",
           "exen": "Exenatide", "dula": "Dulaglutide", "lixi": "Lixisenatide",
           "setm": "Setmelanotide"}
DRUG_COLORS = {"sema": "#c2185b", "lira": "#e65100", "tirz": "#6a1b9a",
               "exen": "#2e7d32", "dula": "#1565c0", "lixi": "#8d6e63",
               "setm": "#00838f"}
EVENT_CN = {"Suicidality": "自杀相关", "Depression": "抑郁", "Anxiety": "焦虑",
            "Eating disorder": "进食障碍", "Psychosis": "精神病性症状"}

byyear = pd.read_csv("results/faers_full_ped_glp1_byyear.csv").fillna(0)
ped = pd.read_csv("results/faers_full_ror_pediatric.csv")
adu = pd.read_csv("results/faers_full_ror_adult.csv")

fig = plt.figure(figsize=(14, 6), dpi=150)
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.35], wspace=0.24)

# ===== Panel A: 年度趋势 =====
ax = fig.add_subplot(gs[0])
years = byyear["yr"].astype(int).tolist()
xl = [str(y) + ("*" if y == 2026 else "") for y in years]
bottom = np.zeros(len(years))
for d in ["sema", "lira", "tirz", "setm", "dula", "exen", "lixi"]:
    v = byyear[d].values.astype(float)
    ax.bar(range(len(years)), v, bottom=bottom, color=DRUG_COLORS[d],
           label=DRUG_CN[d], width=0.72)
    bottom += v
for i, tot in enumerate(byyear["glp_any"].values):
    if tot > 0:
        ax.text(i, bottom[i] + 2, f"{int(tot)}", ha="center", fontsize=8, color="#333333")
ax.set_xticks(range(len(years)))
ax.set_xticklabels(xl, fontsize=9, rotation=45)
ax.set_ylabel("儿科 GLP-1 类报告数 (例)", fontsize=10.5)
ax.set_ylim(0, max(bottom) * 1.18)
ax.legend(frameon=False, fontsize=8.5, loc="upper left")
ax.set_title("A  FAERS 儿科 GLP-1 类报告年度趋势 (2012Q4-2026Q2)", fontsize=11.5, loc="left", pad=8)
ax.spines[["top", "right"]].set_visible(False)
ax.text(0.0, -0.2, "*2026 年仅含 Q1-Q2；按 caseid 去重后 <18 岁报告；图顶数字为去重病例数，\n堆叠色块为分药计数（一例可含多个 GLP-1 药物，分药合计可略大于去重病例数）",
        transform=ax.transAxes, fontsize=8, color="#555555", va="top")

# ===== Panel B: 类级 ROR 森林图 =====
ax2 = fig.add_subplot(gs[1])
events = list(EVENT_CN.keys())
grp_ped = ped[ped["drug"] == "GLP1-class(all)"].set_index("event_group")
grp_adu = adu[adu["drug"] == "GLP1-class(all)"].set_index("event_group")
offsets = [0.16, -0.16]
for gi, ev in enumerate(events):
    base = -gi
    for si, (grp, col, mk, lab) in enumerate([(grp_ped, "#c2185b", "o", "儿科 (<18 岁, N=616,860)"),
                                              (grp_adu, "#1565c0", "s", "成人 (≥18 岁, N=9,457,685)")]):
        r = grp.loc[ev]
        yy = base + offsets[si]
        ax2.errorbar(r["ROR"], yy, xerr=[[r["ROR"] - r["ROR_lo"]], [r["ROR_hi"] - r["ROR"]]],
                     fmt=mk, color=col, ecolor=col, elinewidth=1.4, capsize=3.5,
                     markersize=6.5, label=lab if gi == 0 else None)
        star = " *" if r["ROR_lo"] > 1 or r["ROR_hi"] < 1 else ""
        ax2.text(38, yy, f"a={int(r['a'])}  {r['ROR']:.2f} ({r['ROR_lo']:.2f}-{r['ROR_hi']:.2f}){star}",
                 va="center", fontsize=8.3, color=col)
ax2.axvline(1.0, color="#888888", ls="--", lw=1.0, zorder=0)
ax2.set_xscale("log")
ticks = [0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20]
ax2.set_xticks(ticks)
ax2.get_xaxis().set_major_formatter(matplotlib.ticker.FixedFormatter([str(t) for t in ticks]))
ax2.set_xlim(0.05, 30)
ax2.set_yticks([-i for i in range(len(events))])
ax2.set_yticklabels([EVENT_CN[e] for e in events], fontsize=10.5)
ax2.set_ylim(-len(events) + 0.5, 1.1)
ax2.set_xlabel("ROR (95% CI)，对数刻度；* = 95% CI 不含 1", fontsize=10.5)
ax2.set_title("B  GLP-1 类 × 精神科事件组: 儿科 vs 成人 ROR", fontsize=11.5, loc="left", pad=8)
ax2.legend(frameon=False, fontsize=9, loc="upper left", borderaxespad=0.1)
ax2.spines[["top", "right"]].set_visible(False)
ax2.text(38, 0.75, "a    ROR (95% CI)", fontsize=8.5, color="#333333", va="center")

fig.suptitle("FAERS 全量分析 (2012Q4-2026Q2, 17,585,762 例去重病例): GLP-1 类精神科不相称性",
             fontsize=12.5, y=0.995)
fig.savefig("results/faers_full_forest.png", bbox_inches="tight")
print("saved results/faers_full_forest.png")
