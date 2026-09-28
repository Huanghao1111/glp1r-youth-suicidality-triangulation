#!/usr/bin/env python3
"""Assemble phase-2 zip-pass flags onto deduplicated cases and compute:
1. Primary-suspect-only GLP-1 RA disproportionality (pediatric + adult, all event groups).
2. ADHD comparator class: pediatric suicidality ROR overall and pre/post-2023.
3. GI control events (nausea/vomiting/diarrhoea/constipation/pancreatitis):
   pediatric GLP-1 RA ROR pre/post-2023, contrasted with suicidality.
"""
import math
import duckdb
import pandas as pd
from scipy.stats import gamma as gamma_dist

con = duckdb.connect()
con.execute("pragma memory_limit='8GB'")
con.execute("create or replace view raw as select * from read_parquet('data/faers_proc/*.parquet')")
con.execute("""
create or replace view del as
select distinct try_cast(column0 as bigint) cid
from read_csv('data/faers_proc/deleted_caseids.txt', header=false, all_varchar=true)
where try_cast(column0 as bigint) is not null
""")
con.execute("""
create or replace view drugflags as
select try_cast(pid as bigint) pid, max(glp_ps) glp_ps, max(adhd) adhd
from read_csv('data/faers_flags/*_drug.csv', header=true)
group by pid
""")
con.execute("""
create or replace view giflags as
select distinct try_cast(pid as bigint) pid
from read_csv('data/faers_flags/*_gi.csv', header=true)
""")
con.execute("""
create or replace view cases as
with ranked as (
  select r.*, row_number() over (
    partition by caseid
    order by coalesce(fda_dt,0) desc, coalesce(caseversion,0) desc, primaryid desc
  ) rn
  from raw r
  where caseid not in (select cid from del)
)
select k.primaryid, k.caseid, k.fda_dt, k.quarter, k.age_years,
       k.sema, k.lira, k.tirz, k.exen, k.dula, k.lixi, k.setm,
       case when k.sema+k.lira+k.tirz+k.exen+k.dula+k.lixi+k.setm>0 then 1 else 0 end as glp,
       k.ev_suic, k.ev_depr, k.ev_anx, k.ev_eat, k.ev_psy,
       coalesce(df.glp_ps,0) as glp_ps, coalesce(df.adhd,0) as adhd,
       case when gi.pid is not null then 1 else 0 end as gi
from (select * from ranked where rn=1) k
left join drugflags df on k.primaryid = df.pid
left join giflags gi on k.primaryid = gi.pid
""")

def ror_ci(a, b, c, d):
    corr = 0.5 if min(a, b, c, d) == 0 else 0.0
    a_, b_, c_, d_ = a + corr, b + corr, c + corr, d + corr
    v = (a_ * d_) / (b_ * c_)
    se = math.sqrt(1/a_ + 1/b_ + 1/c_ + 1/d_)
    return v, math.exp(math.log(v) - 1.96*se), math.exp(math.log(v) + 1.96*se), corr > 0

def ic025(a, n_drug, n_event, n_total):
    E = n_drug * n_event / n_total
    lo = math.log2(gamma_dist.ppf(0.025, a + 0.5, scale=1/(1 + E)))
    pt = math.log2((a + 0.5) / (E + 0.5))
    return pt, lo

def cell(drug_cond, event, where):
    q = f"""select count(*) n_total,
              sum(case when {drug_cond} then 1 else 0 end) n_drug,
              sum({event}) n_event,
              sum(case when {drug_cond} then {event} else 0 end) a
            from cases where {where}"""
    n_total, n_drug, n_event, a = (int(x) for x in con.execute(q).fetchone())
    b, c = n_drug - a, n_event - a
    d = n_total - n_drug - n_event + a
    v, lo, hi, corr = ror_ci(a, b, c, d)
    icpt, iclo = ic025(a, n_drug, n_event, n_total)
    return dict(N=n_total, n_drug=n_drug, a=a, ROR=round(v, 3), lo=round(lo, 3),
                hi=round(hi, 3), haldane=corr, IC=round(icpt, 3), IC025=round(iclo, 3))

EVENTS = {"ev_suic": "Suicidality", "ev_depr": "Depression", "ev_anx": "Anxiety",
          "ev_eat": "Eating disorder", "ev_psy": "Psychosis"}

# ---- 1. primary-suspect only ----
rows = []
for stratum, where in [("pediatric", "age_years<18"), ("adult", "age_years>=18")]:
    for ek, el in EVENTS.items():
        r = cell("glp_ps=1", ek, where)
        rows.append(dict(stratum=stratum, exposure="GLP-1 RA, primary suspect only",
                         event=el, **r))
ps = pd.DataFrame(rows)
ps.to_csv("results/faers_ps_only_ror.csv", index=False)
print("=== 1. primary-suspect-only class-level ROR ===")
print(ps[["stratum", "event", "n_drug", "a", "ROR", "lo", "hi", "IC025"]].to_string(index=False))

# ---- 2. ADHD comparator, pediatric suicidality ----
rows = []
for period, where in [("overall", "age_years<18"),
                      ("pre-2023", "age_years<18 and fda_dt<20230101"),
                      ("post-2023", "age_years<18 and fda_dt>=20230101")]:
    r = cell("adhd=1", "ev_suic", where)
    rows.append(dict(period=period, exposure="ADHD comparator class", event="Suicidality", **r))
    r = cell("glp=1", "ev_suic", where)
    rows.append(dict(period=period, exposure="GLP-1 RA (reference)", event="Suicidality", **r))
cmp_df = pd.DataFrame(rows)
cmp_df.to_csv("results/faers_adhd_comparator.csv", index=False)
print("\n=== 2. ADHD comparator vs GLP-1, pediatric suicidality ===")
print(cmp_df[["period", "exposure", "n_drug", "a", "ROR", "lo", "hi", "IC025"]].to_string(index=False))

# ---- 3. GI control events pre/post-2023 ----
rows = []
for period, where in [("pre-2023", "age_years<18 and fda_dt<20230101"),
                      ("post-2023", "age_years<18 and fda_dt>=20230101")]:
    r = cell("glp=1", "gi", where)
    rows.append(dict(period=period, event="GI control (nausea/vomiting/diarrhoea/constipation/pancreatitis)", **r))
    r = cell("glp=1", "ev_suic", where)
    rows.append(dict(period=period, event="Suicidality (reference)", **r))
gi = pd.DataFrame(rows)
gi.to_csv("results/faers_gi_control.csv", index=False)
print("\n=== 3. GI control vs suicidality, pediatric GLP-1 RA ===")
print(gi[["period", "event", "n_drug", "a", "ROR", "lo", "hi", "IC025"]].to_string(index=False))
con.close()
