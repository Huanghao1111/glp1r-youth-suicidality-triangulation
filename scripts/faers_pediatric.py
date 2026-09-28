#!/usr/bin/env python3
"""FAERS pediatric disproportionality analysis via openFDA API. Run: python faers_pediatric.py pediatric|adult"""
import urllib.request, urllib.parse, json, time, math, sys, csv, os

BASE = "https://api.fda.gov/drug/event.json"

DRUGS = {
    "Semaglutide":  "SEMAGLUTIDE",
    "Liraglutide":  "LIRAGLUTIDE",
    "Tirzepatide":  "TIRZEPATIDE",
    "Exenatide":    "EXENATIDE",
    "Dulaglutide":  "DULAGLUTIDE",
    "Lixisenatide": "LIXISENATIDE",
    "Setmelanotide":"SETMELANOTIDE",
}
EVENT_GROUPS = {
    "Suicidality": ["SUICIDAL IDEATION","SUICIDE ATTEMPT","COMPLETED SUICIDE",
                     "INTENTIONAL SELF-INJURY","SELF-INJURIOUS IDEATION","SUICIDAL BEHAVIOUR"],
    "Depression":  ["DEPRESSION","DEPRESSED MOOD","MAJOR DEPRESSION","DEPRESSIVE SYMPTOM"],
    "Anxiety":     ["ANXIETY","ANXIETY DISORDER","PANIC ATTACK","PANIC DISORDER"],
    "Eating disorder": ["ANOREXIA NERVOSA","BULIMIA NERVOSA","BINGE EATING DISORDER",
                         "EATING DISORDER","AVOIDANT/RESTRICTIVE FOOD INTAKE DISORDER"],
    "Psychosis":   ["PSYCHOTIC DISORDER","HALLUCINATION","DELUSION","PSYCHOTIC SYMPTOM"],
}
STRATA = {"pediatric": ["3", "4"], "adult": ["5", "6"]}

def q(search):
    url = BASE + "?search=" + urllib.parse.quote(search, safe='()+:"') + "&limit=1"
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=30) as r:
                d = json.loads(r.read().decode())
            return d["meta"]["results"]["total"]
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return 0
            time.sleep(1.5 + attempt * 2)
        except Exception:
            time.sleep(1.5 + attempt * 2)
    print("FAILED:", search, file=sys.stderr, flush=True)
    return None

def drug_q(name):  return f'patient.drug.openfda.substance_name:"{name}"'
def class_q():     return "(" + "+OR+".join(drug_q(v) for v in DRUGS.values()) + ")"
def event_q(pts):  return "(" + "+OR+".join(f'patient.reaction.reactionmeddrapt.exact:"{p}"' for p in pts) + ")"
def age_q(codes):  return "(" + "+OR+".join(f"patient.patientagegroup:{c}" for c in codes) + ")"

def ror(a, b, c, d):
    corr = 0.5 if min(a, b, c, d) == 0 else 0.0
    a_, b_, c_, d_ = a + corr, b + corr, c + corr, d + corr
    v = (a_ * d_) / (b_ * c_)
    se = math.sqrt(1/a_ + 1/b_ + 1/c_ + 1/d_)
    return v, math.exp(math.log(v) - 1.96 * se), math.exp(math.log(v) + 1.96 * se), corr > 0

stratum = sys.argv[1] if len(sys.argv) > 1 else "pediatric"
ages = STRATA[stratum]
A = age_q(ages)

n_total = q(A)
print(f"[{stratum}] total reports: {n_total}", flush=True)

targets = {**DRUGS, "GLP1-class(all)": None}
n_drug, n_both = {}, {}
for drug, sub in targets.items():
    ds = class_q() if sub is None else drug_q(sub)
    n_drug[drug] = q(f"{A}+AND+{ds}")
    print(f"[{stratum}] {drug} total: {n_drug[drug]}", flush=True)

n_event = {grp: q(f"{A}+AND+{event_q(pts)}") for grp, pts in EVENT_GROUPS.items()}
for grp, n in n_event.items():
    print(f"[{stratum}] event {grp}: {n}", flush=True)

out = f"results/faers_ror_{stratum}.csv"
new = not os.path.exists(out)
with open(out, "a", newline="") as f:
    w = csv.writer(f)
    if new:
        w.writerow(["stratum","drug","event_group","a","b","c","d","ror","lo","hi","haldane","signal"])
    for drug, sub in targets.items():
        ds = class_q() if sub is None else drug_q(sub)
        for grp, pts in EVENT_GROUPS.items():
            a = q(f"{A}+AND+{ds}+AND+{event_q(pts)}")
            if a is None or n_drug[drug] is None or n_event[grp] is None:
                continue
            b = n_drug[drug] - a
            c = n_event[grp] - a
            d = n_total - n_drug[drug] - n_event[grp] + a
            v, lo, hi, corr = ror(a, b, c, d)
            w.writerow([stratum, drug, grp, a, b, c, d, round(v,3), round(lo,3), round(hi,3), corr, lo > 1]); f.flush()
            print(f"{stratum:9s} {drug:16s} {grp:15s} a={a:<4d} ROR={v:5.2f} ({lo:.2f}-{hi:.2f}){'*' if lo>1 else ''}", flush=True)
print("saved ->", out)
