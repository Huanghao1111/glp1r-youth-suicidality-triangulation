#!/usr/bin/env python3
"""FAERS phase-2 analyses on the processed parquet corpus:
A. Quarter-resolved pediatric class-level suicidality ROR/IC025 + report volume.
B. Masking assessment: recompute pediatric suicidality RORs with semaglutide
   removed from the comparator background (drug level) and class level ex-semaglutide.
C. Age-missingness profile by year / GLP-1 / event / seriousness, and sensitivity
   ROR contrasting pediatric reports against ALL non-pediatric reports (incl. unknown age).
Conventions identical to scripts/faers_analyze.py (dedup, Haldane, WHO-UMC IC025).
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
create or replace view cases as
with ranked as (
  select r.*, row_number() over (
    partition by caseid
    order by coalesce(fda_dt,0) desc, coalesce(caseversion,0) desc, primaryid desc
  ) rn
  from raw r
  where caseid not in (select cid from del)
)
select primaryid, caseid, fda_dt, quarter, age_years, sex, occr_country,
       sema, lira, tirz, exen, dula, lixi, setm,
       case when sema+lira+tirz+exen+dula+lixi+setm>0 then 1 else 0 end as glp,
       ev_suic, ev_depr, ev_anx, ev_eat, ev_psy, serious
from ranked where rn=1
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

def cell(drug_where, event, base_where):
    q = f"""select count(*) n_total,
              sum(case when {drug_where} then 1 else 0 end) n_drug,
              sum({event}) n_event,
              sum(case when {drug_where} then {event} else 0 end) a
            from cases where {base_where}"""
    r = con.execute(q).fetchone()
    n_total, n_drug, n_event, a = (int(x) for x in r)
    b, c = n_drug - a, n_event - a
    d = n_total - n_drug - n_event + a
    v, lo, hi, corr = ror_ci(a, b, c, d)
    icpt, iclo = ic025(a, n_drug, n_event, n_total)
    return dict(N=n_total, n_drug=n_drug, n_event=n_event, a=a,
                ROR=round(v, 3), lo=round(lo, 3), hi=round(hi, 3), haldane=corr,
                IC=round(icpt, 3), IC025=round(iclo, 3))

# ---------- A. quarter-resolved pediatric suicidality ----------
qa = con.execute("""
select quarter, count(*) n_ped, sum(glp) n_glp, sum(ev_suic) n_suic,
       sum(glp*ev_suic) a
from cases where age_years<18 group by quarter order by quarter
""").fetchdf()
rows = []
for _, r in qa.iterrows():
    n_total, n_drug, n_event, a = int(r.n_ped), int(r.n_glp), int(r.n_suic), int(r.a)
    b, c = n_drug - a, n_event - a
    d = n_total - n_drug - n_event + a
    v, lo, hi, corr = ror_ci(a, b, c, d)
    icpt, iclo = ic025(a, n_drug, n_event, n_total)
    rows.append([r.quarter, n_total, n_drug, n_event, a,
                 round(v, 2), round(lo, 2), round(hi, 2), corr, round(iclo, 3)])
qt = pd.DataFrame(rows, columns=["quarter", "N_ped", "n_glp", "n_suic", "a",
                                 "ROR", "lo", "hi", "haldane", "IC025"])
qt.to_csv("results/faers_ped_suic_quarterly.csv", index=False)
print("=== A. quarterly pediatric class suicidality (nonzero-a quarters) ===")
print(qt[qt.a > 0].to_string(index=False))

# ---------- B. masking ----------
rows = []
# class level excluding semaglutide from exposure AND background
r = cell("(lira+tirz+exen+dula+lixi+setm)>0", "ev_suic", "age_years<18 and sema=0")
rows.append(dict(analysis="class (ex-semaglutide), background ex-semaglutide", **r))
# class level as reported (for reference)
r = cell("glp=1", "ev_suic", "age_years<18")
rows.append(dict(analysis="class (all), full background [reference]", **r))
# drug level, background excluding semaglutide reports
for d in ["lira", "tirz", "exen", "dula", "lixi", "setm"]:
    r = cell(f"{d}=1", "ev_suic", f"age_years<18 and sema=0")
    rows.append(dict(analysis=f"{d} vs background ex-semaglutide", **r))
# semaglutide vs background excluding ALL GLP-1 (unmasking from class siblings)
r = cell("sema=1", "ev_suic", "age_years<18 and glp=0 or (age_years<18 and sema=1)")
rows.append(dict(analysis="semaglutide vs background ex-all-GLP1", **r))
mk = pd.DataFrame(rows)
mk.to_csv("results/faers_ped_suic_masking.csv", index=False)
print("\n=== B. masking assessment (pediatric suicidality) ===")
print(mk[["analysis", "N", "n_drug", "a", "ROR", "lo", "hi", "IC025"]].to_string(index=False))

# ---------- C. age-missingness profile ----------
prof_year = con.execute("""
select floor(fda_dt/10000)::int yr, count(*) n,
       round(100.0*count(*) filter (age_years is null)/count(*), 1) pct_missing
from cases where fda_dt is not null group by 1 order by 1
""").fetchdf()
prof_glp = con.execute("""
select case when glp=1 then 'GLP-1 RA' else 'all other drugs' end grp,
       count(*) n, round(100.0*count(*) filter (age_years is null)/count(*), 1) pct_missing
from cases group by 1
""").fetchdf()
prof_ev = con.execute("""
select case when ev_suic=1 then 'suicidality report' else 'no suicidality' end grp,
       count(*) n, round(100.0*count(*) filter (age_years is null)/count(*), 1) pct_missing
from cases group by 1
""").fetchdf()
prof_ser = con.execute("""
select case when serious=1 then 'serious' else 'non-serious/unknown' end grp,
       count(*) n, round(100.0*count(*) filter (age_years is null)/count(*), 1) pct_missing
from cases group by 1
""").fetchdf()
prof_glp_year = con.execute("""
select floor(fda_dt/10000)::int yr, glp,
       count(*) n, round(100.0*count(*) filter (age_years is null)/count(*), 1) pct_missing
from cases where fda_dt is not null group by 1, 2 order by 1, 2
""").fetchdf()
prof_year.to_csv("results/faers_agemiss_by_year.csv", index=False)
prof_glp_year.to_csv("results/faers_agemiss_by_year_glp.csv", index=False)
print("\n=== C. age missingness ===")
print("overall by GLP-1:\n", prof_glp.to_string(index=False))
print("by suicidality event:\n", prof_ev.to_string(index=False))
print("by seriousness:\n", prof_ser.to_string(index=False))
print("by year:\n", prof_year.to_string(index=False))

# sensitivity: adult stratum expanded to include age-unknown reports (worst-case comparator)
r1 = cell("glp=1", "ev_suic", "age_years<18")                      # pediatric, primary
r2 = cell("glp=1", "ev_suic", "age_years>=18")                     # adult, primary
r3 = cell("glp=1", "ev_suic", "age_years>=18 or age_years is null")  # adult + unknown-age
r4 = cell("glp=1", "ev_suic", "age_years is null")                 # unknown-age alone
sens = pd.DataFrame([
    dict(analysis="pediatric (<18) [primary]", **r1),
    dict(analysis="adult (>=18) [primary]", **r2),
    dict(analysis="adult + age-unknown [sensitivity]", **r3),
    dict(analysis="age-unknown alone [profile]", **r4),
])
# z contrast pediatric vs expanded adult
import scipy.stats as st
b1, s1 = math.log(r1["ROR"]), (math.log(r1["hi"]) - math.log(r1["lo"])) / 3.92
b3, s3 = math.log(r3["ROR"]), (math.log(r3["hi"]) - math.log(r3["lo"])) / 3.92
z = (b1 - b3) / math.sqrt(s1**2 + s3**2)
print(f"pediatric vs adult+unknown z={z:.2f}, P={2*st.norm.sf(abs(z)):.2e}")
sens.to_csv("results/faers_agemiss_sensitivity.csv", index=False)
print("\n=== C. age-missingness sensitivity ===")
print(sens[["analysis", "N", "n_drug", "a", "ROR", "lo", "hi", "IC025"]].to_string(index=False))
con.close()
