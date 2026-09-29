#!/usr/bin/env python3
"""Phase-3 comparator panel: pediatric suicidality disproportionality pre/post-2023
for SSRI, second-generation antipsychotics, montelukast, isotretinoin comparator classes,
with GLP-1 RA and ADHD class as references. Conventions as in faers_phase2_assemble.py.
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
create or replace view cmpflags as
select try_cast(pid as bigint) pid, max(ssri) ssri, max(sga) sga, max(mont) mont, max(iso) iso
from read_csv('data/faers_flags3/*_cmp.csv', header=true)
group by pid
""")
con.execute("""
create or replace view adhdflags as
select try_cast(pid as bigint) pid, max(adhd) adhd
from read_csv('data/faers_flags/*_drug.csv', header=true)
group by pid
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
select k.primaryid, k.fda_dt, k.age_years,
       case when k.sema+k.lira+k.tirz+k.exen+k.dula+k.lixi+k.setm>0 then 1 else 0 end as glp,
       k.ev_suic,
       coalesce(c.ssri,0) as ssri, coalesce(c.sga,0) as sga,
       coalesce(c.mont,0) as mont, coalesce(c.iso,0) as iso,
       coalesce(af.adhd,0) as adhd
from (select * from ranked where rn=1) k
left join cmpflags c on k.primaryid = c.pid
left join adhdflags af on k.primaryid = af.pid
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
    return lo

def cell(drug_cond, where):
    q = f"""select count(*) n_total,
              sum(case when {drug_cond} then 1 else 0 end) n_drug,
              sum(ev_suic) n_event,
              sum(case when {drug_cond} then ev_suic else 0 end) a
            from cases where {where}"""
    n_total, n_drug, n_event, a = (int(x) for x in con.execute(q).fetchone())
    b, c = n_drug - a, n_event - a
    d = n_total - n_drug - n_event + a
    v, lo, hi, corr = ror_ci(a, b, c, d)
    return dict(N=n_total, n_drug=n_drug, a=a, ROR=round(v, 2), lo=round(lo, 2),
                hi=round(hi, 2), haldane=corr, IC025=round(ic025(a, n_drug, n_event, n_total), 2))

CLASSES = [("glp=1", "GLP-1 RA (index)"),
           ("adhd=1", "ADHD drugs"),
           ("ssri=1", "SSRIs"),
           ("sga=1", "Second-generation antipsychotics"),
           ("mont=1", "Montelukast"),
           ("iso=1", "Isotretinoin")]

rows = []
for cond, lab in CLASSES:
    for period, where in [("pre-2023", "age_years<18 and fda_dt<20230101"),
                          ("post-2023", "age_years<18 and fda_dt>=20230101")]:
        r = cell(cond, where)
        rows.append(dict(drug_class=lab, period=period, **r))
df = pd.DataFrame(rows)
# fold change post/pre per class
piv = df.pivot_table(index="drug_class", columns="period", values="ROR")
piv["fold"] = (piv["post-2023"] / piv["pre-2023"]).round(2)
df.to_csv("results/faers_comparator_panel.csv", index=False)
print(df.to_string(index=False))
print("\n=== fold-change (post/pre) ===")
print(piv[["pre-2023", "post-2023", "fold"]].to_string())
con.close()
