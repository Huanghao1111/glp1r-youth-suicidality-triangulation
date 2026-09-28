#!/usr/bin/env python3
"""YRBS 2025 sensitivity: extend models with available proximal covariates.
2025 questionnaire does NOT include unhealthy weight-control items (fasting/diet pills/
vomiting-laxatives) — documented absence. Available: QN23 school bullying, QN24 electronic
bullying, QN25 non-suicidal self-injury, QN94 binge-eating with loss of control,
Q62 sexual identity, QN85/QN86 food insecurity (SES proxy).
"""
import re
import numpy as np
import pandas as pd
import statsmodels.api as sm

pat = re.compile(r"^@(\d+)\s+(\S+)\s+(\$?)(\d+)(?:\.(\d+))?", re.M)
txt = open("data/2025XXH-SAS-Input-Program.sas", encoding="utf-8", errors="ignore").read()
spec = {m.group(2): (int(m.group(1))-1, int(m.group(1))-1+int(m.group(4)), m.group(3)=="$")
        for m in pat.finditer(txt)}

NEED = ["Q1", "Q2", "QN26", "QN27", "QN29", "Q65", "Q62", "Q64",
        "QN23", "QN24", "QN25", "QN94", "QN85", "QN86",
        "WEIGHT", "STRATUM", "PSU", "BMIPCT", "RACEETH"]
rows = []
with open("data/XXH2025_YRBS_Data.dat", encoding="ascii", errors="ignore") as f:
    for line in f:
        if len(line) < 425:
            continue
        rows.append({v: line[spec[v][0]:spec[v][1]].strip() for v in NEED})
df = pd.DataFrame(rows)
for v in NEED:
    if v != "RACEETH":
        df[v] = pd.to_numeric(df[v], errors="coerce")
df["sex"] = df["Q2"].map({1: "Female", 2: "Male"})
df["raceeth"] = df["RACEETH"].map({"1": "AI/AN", "2": "Asian", "3": "Black", "4": "Hispanic/Latino",
                                   "5": "Mideast/N African", "6": "NH/PI", "7": "White",
                                   "8": "Multiple-Hisp", "9": "Multiple-NonHisp"})
pct = df["BMIPCT"].copy()
if pct.dropna().max() <= 1.0:
    pct *= 100
df["bmi_grp"] = np.where(pct >= 95, "Obese", np.where(pct >= 85, "Overweight",
                 np.where(pct >= 5, "Normal", "Underweight")))
df["trylose"] = (df["Q65"] == 1).where(df["Q65"].isin([1, 2, 3, 4]))
df["perceive_ovw"] = (df["Q64"] >= 4).where(df["Q64"].isin([1, 2, 3, 4, 5]))
for q in ["QN23", "QN24", "QN25", "QN94", "QN85", "QN86", "QN26", "QN27", "QN29"]:
    df[q + "_y"] = (df[q] == 1).where(df[q].isin([1, 2]))
df["lgb"] = np.where(df["Q62"].isin([2, 3, 4]), 1, np.where(df["Q62"] == 1, 0, np.nan))
df["id_unsure"] = np.where(df["Q62"].isin([5, 6]), 1, np.where(df["Q62"].isin([1, 2, 3, 4]), 0, np.nan))

ana = df.dropna(subset=["WEIGHT", "STRATUM", "PSU", "sex"]).copy()
ana = ana[ana["bmi_grp"].isin(["Normal", "Overweight", "Obese"])]

BASE = ["sex", "raceeth", "Q1"]
EXT = ["QN23_y", "QN24_y", "QN25_y", "QN94_y", "QN85_y", "lgb", "id_unsure"]

def fit(d, y, exp, covs, label, model):
    d = d.dropna(subset=[y, exp] + covs).reset_index(drop=True)
    if d[exp].nunique() < 2:
        return None
    parts = []
    if "sex" in covs:
        parts.append(pd.get_dummies(d[["sex", "raceeth"]], drop_first=True))
    X = pd.DataFrame({"exp": d[exp].astype(float)})
    if "Q1" in covs:
        X["age"] = d["Q1"].astype(float)
    if parts:
        X = pd.concat([X, parts[0]], axis=1)
    for c in covs:
        if c not in ("sex", "raceeth", "Q1"):
            X[c] = d[c].astype(float)
    X = sm.add_constant(X.astype(float))
    nz = [c for c in X.columns if c == "const" or X[c].std() > 0]
    X = X[nz]
    res = sm.GLM(d[y].astype(float), X, family=sm.families.Binomial(),
                 freq_weights=d["WEIGHT"]).fit(cov_type="cluster",
                                               cov_kwds={"groups": d["PSU"].values})
    b, se, p = res.params["exp"], res.bse["exp"], res.pvalues["exp"]
    return dict(analysis=label, model=model, n=len(d),
                OR=np.exp(b), lo=np.exp(b-1.96*se), hi=np.exp(b+1.96*se), p=p)

OUT = {"QN26_y": "Persistent sadness/hopelessness",
       "QN27_y": "Seriously considered suicide",
       "QN29_y": "Suicide attempt"}
rows_out = []
for y, lab in OUT.items():
    for exp, elab in [("trylose", "Trying to lose weight"), ("perceive_ovw", "Perceived overweight (normal BMI)")]:
        dd = ana if exp == "trylose" else ana[ana.bmi_grp == "Normal"]
        r1 = fit(dd, y, exp, BASE, f"{elab} -> {lab}", "Base (age, sex, race/ethnicity)")
        r2 = fit(dd, y, exp, BASE + EXT, f"{elab} -> {lab}",
                 "Extended (+ school bullying, electronic bullying, NSSI, binge-eating LOC, food insecurity, sexual identity)")
        for r in (r1, r2):
            if r:
                rows_out.append(r)
    # stratum-restricted normal-weight trylose, extended only vs base
    dd = ana[ana.bmi_grp == "Normal"]
    r1 = fit(dd, y, "trylose", BASE, f"Trying to lose weight -> {lab} (normal-weight only)", "Base (age, sex, race/ethnicity)")
    r2 = fit(dd, y, "trylose", BASE + EXT, f"Trying to lose weight -> {lab} (normal-weight only)",
             "Extended (+ school bullying, electronic bullying, NSSI, binge-eating LOC, food insecurity, sexual identity)")
    for r in (r1, r2):
        if r:
            rows_out.append(r)

res = pd.DataFrame(rows_out)
res.to_csv("results/yrbs_2025_extended_covariates.csv", index=False)
pd.set_option("display.width", 200)
print(res.round(3).to_string(index=False))
