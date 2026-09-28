# -*- coding: utf-8 -*-
# 2025 YRBS 减重行为 × 心理健康: A=减重意向流行率(BMI×性别) B=OR森林图
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.sans-serif"] = ["Noto Sans CJK SC", "SimHei", "Microsoft YaHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

# ---------- Panel A data ----------
prev = {}  # (sex, bmi) -> (pct, lo, hi, n)
with open("results/yrbs_2025_trylose_bysex.csv", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        prev[(r["sex"], r["bmi"])] = (float(r["pct"]), float(r["lo"]), float(r["hi"]), int(r["n"]))

# ---------- Panel B data ----------
rows = []
with open("results/yrbs_2025_weightcontrol_or.csv", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        rows.append((r["analysis"], float(r["OR"]), float(r["lo"]), float(r["hi"]), float(r["p"])))

def pick(sub):
    out = {}
    for name, OR, lo, hi, p in rows:
        if sub in name:
            if "全体" in name:
                out["all"] = (OR, lo, hi)
            elif "限Normal" in name:
                out["normal"] = (OR, lo, hi)
            elif "感知超重" in name:
                out["misp"] = (OR, lo, hi)
    return out

outcomes = ["持续悲伤/绝望", "认真考虑过自杀", "自杀未遂"]
data = {o: pick(o) for o in outcomes}

# ---------- figure ----------
fig = plt.figure(figsize=(13.5, 5.6), dpi=150)
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.25], wspace=0.28)

# ===== Panel A =====
ax = fig.add_subplot(gs[0])
bmis = ["Normal", "Overweight", "Obese"]
bmi_cn = {"Normal": "正常体重", "Overweight": "超重", "Obese": "肥胖"}
sexes = ["Female", "Male"]
colors = {"Female": "#c2185b", "Male": "#1565c0"}
x = np.arange(len(bmis)); w = 0.36
truncated = False
for i, sx in enumerate(sexes):
    vals, errs_lo, errs_hi = [], [], []
    for b in bmis:
        pct, lo, hi, n = prev[(sx, b)]
        hi_c = min(hi, 100.0)
        if hi > 100.0:
            truncated = True
        vals.append(pct)
        errs_lo.append(max(pct - lo, 0))
        errs_hi.append(max(hi_c - pct, 0))
    bars = ax.bar(x + (i - 0.5) * w, vals, w, color=colors[sx], alpha=0.88,
                  label="女生" if sx == "Female" else "男生",
                  yerr=[errs_lo, errs_hi], capsize=3.5,
                  error_kw=dict(lw=1.1, ecolor="#444444"))
    for xi, v in zip(x + (i - 0.5) * w, vals):
        ax.text(xi, v + 6.5, f"{v:.1f}", ha="center", va="bottom", fontsize=8.5, color="#333333")
ax.set_xticks(x)
ax.set_xticklabels([bmi_cn[b] for b in bmis], fontsize=10.5)
ax.set_ylabel("正在尝试减重的比例 (%)", fontsize=10.5)
ax.set_ylim(0, 118)
ax.set_yticks(range(0, 101, 20))
ax.legend(frameon=False, fontsize=10, loc="upper left")
ax.set_title("A  不同 BMI 分组中正在尝试减重的青少年比例", fontsize=11.5, loc="left", pad=8)
ax.spines[["top", "right"]].set_visible(False)
note = "2025 美国全国 YRBS，加权估计；误差线为 95% CI"
if truncated:
    note += "\n（小域 Taylor SE 致个别 CI 上限 >100%，图中截断于 100%）"
ax.text(0.0, -0.16, note, transform=ax.transAxes, fontsize=8, color="#555555", va="top")

# ===== Panel B =====
ax2 = fig.add_subplot(gs[1])
series = [("all",   "全体（校正年龄/性别/种族）", "#2e7d32", "o"),
          ("normal", "仅限正常体重者",            "#e65100", "s"),
          ("misp",  "感知错位：体重正常但自认超重", "#6a1b9a", "D")]
ypos = {}
y = 0.0
ylabels, yticks = [], []
off = {-1: 0.22, 0: 0.0, 1: -0.22}
for gi, o in enumerate(outcomes):
    base = -gi * 1.6
    ylabels.append(o); yticks.append(base)
    for si, (key, lab, col, mk) in enumerate(series):
        OR, lo, hi = data[o][key]
        yy = base + off[si - 1] if False else base + [0.22, 0.0, -0.22][si]
        ax2.errorbar(OR, yy, xerr=[[OR - lo], [hi - OR]], fmt=mk, color=col,
                     ecolor=col, elinewidth=1.4, capsize=3.5, markersize=6.5,
                     label=lab if gi == 0 else None)
        ax2.text(4.35, yy, f"{OR:.2f} ({lo:.2f}-{hi:.2f})", va="center",
                 fontsize=8.3, color=col)
    if gi < len(outcomes) - 1:
        ax2.axhline(base - 0.8, color="#dddddd", lw=0.8, zorder=0)
ax2.axvline(1.0, color="#888888", ls="--", lw=1.0, zorder=0)
ax2.set_xscale("log")
ax2.set_xticks([0.8, 1, 1.5, 2, 3, 4])
ax2.get_xaxis().set_major_formatter(matplotlib.ticker.FixedFormatter(["0.8", "1", "1.5", "2", "3", "4"]))
ax2.set_xlim(0.8, 4.3)
ax2.set_yticks(yticks)
ax2.set_yticklabels(ylabels, fontsize=10.5)
ax2.set_ylim(min(yticks) - 0.55, max(yticks) + 1.85)  # 顶部留白放图例
ax2.set_xlabel("校正后 OR (95% CI)，对数刻度", fontsize=10.5)
ax2.set_title("B  减重意向 / 体重感知错位 与心理健康结局的关联", fontsize=11.5, loc="left", pad=8)
ax2.legend(frameon=False, fontsize=8.8, loc="upper left", borderaxespad=0.1)
ax2.spines[["top", "right"]].set_visible(False)
ax2.text(4.35, yticks[0] + 0.55, "OR (95% CI)", fontsize=8.5, color="#333333", va="center")

fig.suptitle("美国青少年减重行为与心理健康：2025 全国 YRBS (n≈9,443，加权 + PSU 聚类稳健 SE)",
             fontsize=12.5, y=0.995)
fig.savefig("results/yrbs_2025_weightcontrol.png", bbox_inches="tight")
print("saved results/yrbs_2025_weightcontrol.png")
