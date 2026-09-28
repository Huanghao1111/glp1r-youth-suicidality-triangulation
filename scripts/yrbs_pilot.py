#!/usr/bin/env python3
"""YRBS 2025 national pilot: obesity/overweight x depressive symptoms & suicidality.
Survey-weighted (WEIGHT/STRATUM/PSU). Data: data/XXH2025_YRBS_Data.dat
Layout parsed from data/2025XXH-SAS-Input-Program.sas
"""
import re, os
import numpy as np
import pandas as pd
import statsmodels.api as sm

# ---- parse SAS input program for column specs ----
spec = {}
pat = re.compile(r"^@(\d+)\s+(\S+)\s+(\$?)(\d+)(?:\.(\d+))?", re.M)
txt = open("data/2025XXH-SAS-Input-Program.sas", encoding="utf-8", errors="ignore").read()
for m in pat.finditer(txt):
    col, var, ch, width = int(m.group(1)), m.group(2), m.group(3), int(m.group(4))
    spec[var] = (col - 1, col - 1 + width, ch == "$")

NEED = ["Q1","Q2","Q3","Q26","Q27","Q28","Q29","Q30","QN26","QN27","QN28","QN29","QN30",
        "Q64","Q65","QNOBESE","QNOWT","WEIGHT","STRATUM","PSU","BMIPCT","RACEETH"]
missing = [v for v in NEED if v not in spec]
assert not missing, f"missing vars: {missing}"

rows = []
with open("data/XXH2025_YRBS_Data.dat", encoding="ascii", errors="ignore") as f:
    for line in f:
        if len(line) < 425:
            continue
        rec = {}
        for v in NEED:
            s, e, ch = spec[v]
            rec[v] = line[s:e].strip()
        rows.append(rec)
df = pd.DataFrame(rows)
print("raw rows:", len(df))

for v in NEED:
    if v == "RACEETH":
        continue
    df[v] = pd.to_numeric(df[v], errors="coerce")

# ---- recode ----
df["sex"] = df["Q2"].map({1: "Female", 2: "Male"})
df["raceeth"] = df["RACEETH"].map({"1":"AI/AN","2":"Asian","3":"Black","4":"Hispanic/Latino",
                                    "5":"Mideast/N African","6":"NH/PI","7":"White",
                                    "8":"Multiple-Hisp","9":"Multiple-NonHisp"})
for q in ["QN26","QN27","QN28","QN29","QN30","QNOBESE","QNOWT"]:
    df[q+"_y"] = (df[q] == 1).where(df[q].isin([1, 2]))
# BMI groups from percentile (verify scale)
pct = df["BMIPCT"].copy()
if pct.dropna().max() <= 1.0:
    pct = pct * 100
df["bmi_grp"] = np.where(pct >= 95, "Obese", np.where(pct >= 85, "Overweight", np.where(pct >= 5, "Normal", "Underweight")))
# sanity check vs QNOBESE
chk = pd.crosstab(df["bmi_grp"], df["QNOBESE"], dropna=False)
print(chk)

OUT = {"QN26": "持续悲伤/绝望", "QN27": "认真考虑过自杀", "QN28": "制定自杀计划",
       "QN29": "自杀未遂", "QN30": "需就医的自杀未遂"}

ana = df.dropna(subset=["WEIGHT","STRATUM","PSU","bmi_grp","sex"]).copy()
ana = ana[ana["bmi_grp"].isin(["Obese","Overweight","Normal"])]
print("analytic n:", len(ana))

# ---- design-based domain prevalence (Taylor linearization, with-replacement) ----
def wt_prop(d, y):
    d = d.dropna(subset=[y])
    if len(d) == 0: return np.nan, np.nan, 0
    w, ys, st, ps = d["WEIGHT"].values, d[y].values, d["STRATUM"].values, d["PSU"].values
    tot = w.sum()
    p = (w * ys).sum() / tot
    var = 0.0
    for s in np.unique(st):
        m = st == s
        psus = np.unique(ps[m])
        zs = np.array([(w[m & (ps == u)] * ys[m & (ps == u)]).sum() for u in psus])
        if len(zs) > 1:
            var += len(zs) / (len(zs) - 1) * ((zs - zs.mean()) ** 2).sum()
    return p, np.sqrt(var) / tot, len(d)

prev_rows = []
for grp in ["Normal","Overweight","Obese"]:
    d = ana[ana["bmi_grp"] == grp]
    for q, lab in OUT.items():
        p, se, n = wt_prop(d, q + "_y")
        prev_rows.append(dict(bmi=grp, outcome=lab, n=n, pct=100*p, se=100*se,
                              lo=100*(p-1.96*se), hi=100*(p+1.96*se)))
prev = pd.DataFrame(prev_rows)
os.makedirs("results", exist_ok=True)
prev.to_csv("results/yrbs_2025_prevalence.csv", index=False)
print(prev.round(2).to_string(index=False))

# ---- weighted logistic regression (cluster-robust by PSU), adjusted ----
def fit_or(d, y, label):
    d = d.dropna(subset=[y, "sex", "raceeth", "Q1"])
    X = pd.get_dummies(d[["bmi_grp","sex","raceeth"]], drop_first=True)
    X = X.rename(columns={"bmi_grp_Overweight":"Overweight","bmi_grp_Obese":"Obese"})
    keep = ["Overweight","Obese"]
    X["age"] = d["Q1"]
    X = sm.add_constant(X[keep + ["age"] + [c for c in X.columns if c.startswith(("sex_","raceeth_"))]].astype(float))
    mod = sm.GLM(d[y].astype(float), X, family=sm.families.Binomial(),
                 freq_weights=d["WEIGHT"])
    res = mod.fit(cov_type="cluster", cov_kwds={"groups": d["PSU"].values})
    out = []
    for k in keep:
        b, se = res.params[k], res.bse[k]
        out.append(dict(outcome=label, term=k, OR=np.exp(b),
                        lo=np.exp(b-1.96*se), hi=np.exp(b+1.96*se), p=res.pvalues[k]))
    return out

or_rows = []
for q, lab in OUT.items():
    or_rows += fit_or(ana, q + "_y", lab)
orr = pd.DataFrame(or_rows)
orr.to_csv("results/yrbs_2025_or.csv", index=False)
print("\nAdjusted ORs (ref=Normal weight, adj age/sex/raceeth, cluster-robust by PSU):")
print(orr.round(3).to_string(index=False))

# ---- by sex for the three primary outcomes ----
sex_rows = []
for sx in ["Female","Male"]:
    d = ana[ana["sex"] == sx]
    for q in ["QN26","QN27","QN29"]:
        for grp in ["Normal","Overweight","Obese"]:
            dd = d[d["bmi_grp"] == grp]
            p, se, n = wt_prop(dd, q + "_y")
            sex_rows.append(dict(sex=sx, bmi=grp, outcome=OUT[q], n=n, pct=100*p, lo=100*(p-1.96*se), hi=100*(p+1.96*se)))
pd.DataFrame(sex_rows).to_csv("results/yrbs_2025_bysex.csv", index=False)
print("\nsaved: results/yrbs_2025_prevalence.csv, yrbs_2025_or.csv, yrbs_2025_bysex.csv")
