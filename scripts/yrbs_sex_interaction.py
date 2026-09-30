#!/usr/bin/env python3
"""YRBS 2025: sex differences in the weight-loss-intention -> suicidality association.
1) Sex-stratified survey-weighted adjusted ORs (age + race/ethnicity), overall and within BMI strata.
2) Formal sex x trying-to-lose-weight interaction test in the full analytic sample.
Same data pipeline and model family as yrbs_weightcontrol.py."""
import re
import numpy as np
import pandas as pd
import statsmodels.api as sm

pat = re.compile(r"^@(\d+)\s+(\S+)\s+(\$?)(\d+)(?:\.(\d+))?", re.M)
txt = open("data/2025XXH-SAS-Input-Program.sas", encoding="utf-8", errors="ignore").read()
spec = {m.group(2): (int(m.group(1))-1, int(m.group(1))-1+int(m.group(4)), m.group(3)=="$")
        for m in pat.finditer(txt)}

NEED = ["Q1","Q2","QN26","QN27","QN29","Q64","Q65",
        "WEIGHT","STRATUM","PSU","BMIPCT","RACEETH"]
rows = []
with open("data/XXH2025_YRBS_Data.dat", encoding="ascii", errors="ignore") as f:
    for line in f:
        if len(line) < 425: continue
        rows.append({v: line[spec[v][0]:spec[v][1]].strip() for v in NEED})
df = pd.DataFrame(rows)
for v in NEED:
    if v != "RACEETH": df[v] = pd.to_numeric(df[v], errors="coerce")

df["sex"] = df["Q2"].map({1:"Female",2:"Male"})
df["raceeth"] = df["RACEETH"].map({"1":"AI/AN","2":"Asian","3":"Black","4":"Hispanic/Latino",
                                    "5":"Mideast/N African","6":"NH/PI","7":"White",
                                    "8":"Multiple-Hisp","9":"Multiple-NonHisp"})
for q in ["QN26","QN27","QN29"]:
    df[q+"_y"] = (df[q]==1).where(df[q].isin([1,2]))
pct = df["BMIPCT"].copy()
if pct.dropna().max() <= 1.0: pct *= 100
df["bmi_grp"] = np.where(pct>=95,"Obese",np.where(pct>=85,"Overweight",np.where(pct>=5,"Normal","Underweight")))
df["trylose"] = (df["Q65"]==1).where(df["Q65"].isin([1,2,3,4]))
df["perceive_ovw"] = (df["Q64"]>=4).where(df["Q64"].isin([1,2,3,4,5]))

ana = df.dropna(subset=["WEIGHT","STRATUM","PSU","sex"]).copy()
ana = ana[ana["bmi_grp"].isin(["Normal","Overweight","Obese"])]
OUT = {"QN26":"persistent_sadness","QN27":"seriously_considered","QN29":"suicide_attempt"}

def design(d, y, exposure, with_sex=True, interaction=False):
    X = pd.get_dummies(d[["raceeth"]], drop_first=True)
    cols = [X]
    base = pd.DataFrame({"exp": d[exposure].astype(float), "age": d["Q1"].astype(float)}, index=d.index)
    if with_sex:
        sx = pd.get_dummies(d[["sex"]], drop_first=True).astype(float)
        base = pd.concat([base, sx], axis=1)
        if interaction:
            base["exp_x_sex"] = d[exposure].astype(float) * d["sex"].map({"Female":0,"Male":1}).astype(float)
    X = pd.concat([base, X.astype(float)], axis=1)
    return sm.add_constant(X.astype(float))

def fit(d, y, exposure, with_sex=True, interaction=False):
    need = [y, exposure, "raceeth", "Q1"] + (["sex"] if with_sex else [])
    d = d.dropna(subset=need)
    if d[exposure].nunique() < 2: return None
    X = design(d, y, exposure, with_sex, interaction)
    res = sm.GLM(d[y].astype(float), X, family=sm.families.Binomial(), freq_weights=d["WEIGHT"])\
            .fit(cov_type="cluster", cov_kwds={"groups": d["PSU"].values})
    return res, len(d)

rows_out = []
# --- A. sex-stratified ORs, overall + within Normal BMI (headline stratum) + other strata
for sx in ["Female","Male"]:
    sub = ana[ana.sex==sx]
    for q, lab in OUT.items():
        r = fit(sub, q+"_y", "trylose", with_sex=False)
        if r:
            res, n = r; b, se, p = res.params["exp"], res.bse["exp"], res.pvalues["exp"]
            rows_out.append(dict(section="stratified", sex=sx, stratum="all_BMI", outcome=lab,
                                 n=n, OR=np.exp(b), lo=np.exp(b-1.96*se), hi=np.exp(b+1.96*se), p=p))
        for g in ["Normal","Overweight","Obese"]:
            r = fit(sub[sub.bmi_grp==g], q+"_y", "trylose", with_sex=False)
            if r:
                res, n = r; b, se, p = res.params["exp"], res.bse["exp"], res.pvalues["exp"]
                rows_out.append(dict(section="stratified", sex=sx, stratum=g, outcome=lab,
                                     n=n, OR=np.exp(b), lo=np.exp(b-1.96*se), hi=np.exp(b+1.96*se), p=p))
        # discordance among normal-weight, by sex
        dn = sub[sub.bmi_grp=="Normal"]
        r = fit(dn, q+"_y", "perceive_ovw", with_sex=False)
        if r:
            res, n = r; b, se, p = res.params["exp"], res.bse["exp"], res.pvalues["exp"]
            rows_out.append(dict(section="discordance_normalBMI", sex=sx, stratum="Normal", outcome=lab,
                                 n=n, OR=np.exp(b), lo=np.exp(b-1.96*se), hi=np.exp(b+1.96*se), p=p))

# --- B. formal interaction test in full sample (and within Normal BMI)
for q, lab in OUT.items():
    for scope, d0 in [("all_BMI", ana), ("Normal_only", ana[ana.bmi_grp=="Normal"])]:
        r = fit(d0, q+"_y", "trylose", with_sex=True, interaction=True)
        if r:
            res, n = r
            b, se, p = res.params["exp_x_sex"], res.bse["exp_x_sex"], res.pvalues["exp_x_sex"]
            rows_out.append(dict(section="interaction", sex="Male_vs_Female", stratum=scope, outcome=lab,
                                 n=n, OR=np.exp(b), lo=np.exp(b-1.96*se), hi=np.exp(b+1.96*se), p=p))
    # discordance (perceived overweight) x sex interaction among normal-BMI teens
    dn0 = ana[ana.bmi_grp=="Normal"]
    r = fit(dn0, q+"_y", "perceive_ovw", with_sex=True, interaction=True)
    if r:
        res, n = r
        b, se, p = res.params["exp_x_sex"], res.bse["exp_x_sex"], res.pvalues["exp_x_sex"]
        rows_out.append(dict(section="interaction_discordance", sex="Male_vs_Female", stratum="Normal_only", outcome=lab,
                             n=n, OR=np.exp(b), lo=np.exp(b-1.96*se), hi=np.exp(b+1.96*se), p=p))

out = pd.DataFrame(rows_out)
out.to_csv("results/yrbs_2025_sex_interaction.csv", index=False)
pd.set_option("display.width", 200)
print(out.round(3).to_string(index=False))
