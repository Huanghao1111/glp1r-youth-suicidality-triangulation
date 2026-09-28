#!/usr/bin/env python3
"""YRBS 2025: weight perception / trying-to-lose-weight x mental health outcomes.
Reuses layout parsing from yrbs_pilot. Survey-weighted; PSU-cluster robust GLM."""
import re, os
import numpy as np
import pandas as pd
import statsmodels.api as sm

pat = re.compile(r"^@(\d+)\s+(\S+)\s+(\$?)(\d+)(?:\.(\d+))?", re.M)
txt = open("data/2025XXH-SAS-Input-Program.sas", encoding="utf-8", errors="ignore").read()
spec = {m.group(2): (int(m.group(1))-1, int(m.group(1))-1+int(m.group(4)), m.group(3)=="$")
        for m in pat.finditer(txt)}

NEED = ["Q1","Q2","Q26","QN26","QN27","QN28","QN29","QN30","Q64","Q65","QN65",
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
for q in ["QN26","QN27","QN28","QN29","QN30"]:
    df[q+"_y"] = (df[q]==1).where(df[q].isin([1,2]))
pct = df["BMIPCT"].copy()
if pct.dropna().max() <= 1.0: pct *= 100
df["bmi_grp"] = np.where(pct>=95,"Obese",np.where(pct>=85,"Overweight",np.where(pct>=5,"Normal","Underweight")))
df["trylose"] = (df["Q65"]==1).where(df["Q65"].isin([1,2,3,4]))
df["perceive_ovw"] = (df["Q64"]>=4).where(df["Q64"].isin([1,2,3,4,5]))

ana = df.dropna(subset=["WEIGHT","STRATUM","PSU","sex"]).copy()
ana = ana[ana["bmi_grp"].isin(["Normal","Overweight","Obese"])]
OUT = {"QN26":"持续悲伤/绝望","QN27":"认真考虑过自杀","QN29":"自杀未遂"}

def wt_prop(d, y):
    d = d.dropna(subset=[y])
    if len(d)==0: return np.nan, np.nan, 0
    w, ys, st, ps = d["WEIGHT"].values, d[y].values, d["STRATUM"].values, d["PSU"].values
    p = (w*ys).sum()/w.sum()
    var = 0.0
    for s in np.unique(st):
        m = st==s; psus = np.unique(ps[m])
        zs = np.array([(w[m&(ps==u)]*ys[m&(ps==u)]).sum() for u in psus])
        if len(zs)>1: var += len(zs)/(len(zs)-1)*((zs-zs.mean())**2).sum()
    return p, np.sqrt(var)/w.sum(), len(d)

def fit_or(d, y, exposure, label):
    d = d.dropna(subset=[y, exposure, "sex","raceeth","Q1"])
    if d[exposure].nunique() < 2: return []
    X = pd.get_dummies(d[["sex","raceeth"]], drop_first=True)
    X["exp"] = d[exposure].astype(float); X["age"] = d["Q1"].astype(float)
    X = sm.add_constant(pd.concat([X[["exp","age"]], X.filter(regex="^(sex_|raceeth_)")], axis=1).astype(float))
    res = sm.GLM(d[y].astype(float), X, family=sm.families.Binomial(), freq_weights=d["WEIGHT"])\
            .fit(cov_type="cluster", cov_kwds={"groups": d["PSU"].values})
    b, se, p = res.params["exp"], res.bse["exp"], res.pvalues["exp"]
    return [dict(analysis=label, OR=np.exp(b), lo=np.exp(b-1.96*se), hi=np.exp(b+1.96*se), p=p)]

# 1. prevalence: trying to lose weight by BMI x sex
r1 = []
for sx in ["Female","Male"]:
    for g in ["Normal","Overweight","Obese"]:
        p, se, n = wt_prop(ana[(ana.sex==sx)&(ana.bmi_grp==g)], "trylose")
        r1.append(dict(sex=sx, bmi=g, n=n, pct=100*p, lo=100*(p-1.96*se), hi=100*(p+1.96*se)))
prev_tl = pd.DataFrame(r1); prev_tl.to_csv("results/yrbs_2025_trylose_bysex.csv", index=False)
print("== trying to lose weight (weighted %) =="); print(prev_tl.round(1).to_string(index=False))

# 2. discordance among Normal BMI
dn = ana[ana.bmi_grp=="Normal"]
p, se, n = wt_prop(dn, "perceive_ovw")
print(f"\nnormal-BMI teens perceiving self as overweight: {100*p:.1f}% (n={n})")

# 3. ORs: trying-to-lose -> outcomes, overall & within BMI strata; discordance within Normal
orr = []
for q, lab in OUT.items():
    orr += fit_or(ana, q+"_y", "trylose", f"减重意向→{lab} (全体, 校正)")
    for g in ["Normal","Overweight","Obese"]:
        orr += fit_or(ana[ana.bmi_grp==g], q+"_y", "trylose", f"减重意向→{lab} (限{g})")
    orr += fit_or(dn, q+"y" if False else q+"_y", "perceive_ovw", f"感知超重(实际正常)→{lab}")
orr = pd.DataFrame(orr); orr.to_csv("results/yrbs_2025_weightcontrol_or.csv", index=False)
print("\n== adjusted ORs =="); print(orr.round(2).to_string(index=False))
