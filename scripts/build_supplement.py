# -*- coding: utf-8 -*-
"""Build Supplementary Appendix docx (three-line tables, SCI style)."""
import pandas as pd
import numpy as np
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

RES = "results"
MR = "results/mr_cis"

doc = Document()
st = doc.styles["Normal"]
st.font.name = "Arial"; st.font.size = Pt(10)
st.element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")
for sec in doc.sections:
    sec.page_width, sec.page_height = Inches(8.5), Inches(11)

def para(text, bold=False, size=10, space_after=6, italic=False):
    p = doc.add_paragraph()
    r = p.add_run(text); r.bold = bold; r.italic = italic
    r.font.size = Pt(size)
    p.paragraph_format.space_after = Pt(space_after)
    return p

def heading(text, size=12):
    p = para(text, bold=True, size=size, space_after=4)
    return p

_table_counter = {"n": 0}
def caption(label, text):
    _table_counter["n"] += 1
    p = doc.add_paragraph()
    r = p.add_run("Table S%d. " % _table_counter["n"]); r.bold = True; r.font.size = Pt(9)
    r2 = p.add_run(text); r2.font.size = Pt(9)
    p.paragraph_format.space_before = Pt(10); p.paragraph_format.space_after = Pt(3)

def footnote(text):
    p = para(text, size=8, space_after=10)
    for r in p.runs: r.font.color.rgb = RGBColor(0x40, 0x40, 0x40)

def set_border(el_pr, edge, val, sz):
    borders = el_pr.find(qn("w:tblBorders")) if el_pr.tag == qn("w:tblPr") else el_pr.find(qn("w:tcBorders"))
    tag = "w:tblBorders" if el_pr.tag == qn("w:tblPr") else "w:tcBorders"
    if borders is None:
        borders = OxmlElement(tag); el_pr.append(borders)
    e = borders.find(qn("w:" + edge))
    if e is None:
        e = OxmlElement("w:" + edge); borders.append(e)
    e.set(qn("w:val"), val)
    if val != "nil":
        e.set(qn("w:sz"), str(sz)); e.set(qn("w:space"), "0"); e.set(qn("w:color"), "auto")

def three_line(table, header=True):
    tblPr = table._tbl.tblPr
    for edge in ("left", "right", "insideH", "insideV"):
        set_border(tblPr, edge, "nil", 0)
    set_border(tblPr, "top", "single", 12)
    set_border(tblPr, "bottom", "single", 12)
    if header:
        for cell in table.rows[0].cells:
            tcPr = cell._tc.get_or_add_tcPr()
            set_border(tcPr, "bottom", "single", 6)
            for p in cell.paragraphs:
                for r in p.runs: r.bold = True

def add_table(df, widths=None, font=8.5):
    t = doc.add_table(rows=len(df) + 1, cols=len(df.columns))
    t.autofit = True
    for j, c in enumerate(df.columns):
        cell = t.rows[0].cells[j]
        cell.text = ""
        r = cell.paragraphs[0].add_run(str(c)); r.font.size = Pt(font)
    for i in range(len(df)):
        for j, c in enumerate(df.columns):
            cell = t.rows[i + 1].cells[j]
            cell.text = ""
            v = df.iloc[i, j]
            r = cell.paragraphs[0].add_run("" if pd.isna(v) else str(v))
            r.font.size = Pt(font)
    for row in t.rows:
        for cell in row.cells:
            for p in cell.paragraphs:
                p.paragraph_format.space_after = Pt(1)
    three_line(t)
    return t

def f_or(or_, lo, hi, d=2):
    return f"{or_:.{d}f} ({lo:.{d}f}\u2013{hi:.{d}f})"

def f_p(p):
    if pd.isna(p): return ""
    if p < 0.001: return f"{p:.1e}".replace("e-0", "\u00d710\u207b") if False else f"{p:.1e}"
    return f"{p:.3f}".lstrip("0") if p < 1 else f"{p:.3f}"

# ================= title =================
para("Supplementary Appendix", bold=True, size=16, space_after=10)
para("for: GLP-1 Receptor Agonists and Suicide-Attempt Risk in Youth: Triangulation Across "
     "Pharmacovigilance, Adolescent Epidemiology, and Developmental Mendelian Randomization",
     italic=True, size=10, space_after=10)
para("Contents: Tables S1–S3. FAERS pipeline, event-group definitions, and case-level QC. "
     "Tables S4–S8. FAERS full disproportionality results. Tables S9–S13. YRBS 2025 full model tables. "
     "Tables S14–S22. Mendelian randomization supporting analyses. Figure S1. Wave-specific forest plot.",
     size=9, space_after=14)

# ================= S1 =================
heading("S1. FAERS pipeline, event-group definitions, and case-level quality control")

caption("Table S1.1.", "FAERS data processing pipeline (2012Q4\u20132026Q2, 55 quarters).")
flow = pd.read_csv(f"{RES}/faers_full_flow.csv")
lab = {"raw_rows_55q": "Report rows received (55 quarters)",
       "rows_deleted_caseids": "Rows removed per FDA deleted-case list",
       "cases_after_dedup": "Unique cases after deduplication (latest version per case)",
       "pediatric_lt18": "Pediatric stratum (age <18 years)",
       "adult_ge18": "Adult stratum (age \u226518 years)",
       "age_unknown": "Age unknown (excluded from age-stratified analyses)"}
f1 = pd.DataFrame({"Step": [lab.get(i, i) for i in flow.item],
                   "N": [f"{int(v):,}" for v in flow.value]})
add_table(f1)
footnote("Deduplication followed FDA recommendations: the most recent case version per CASEID was retained, "
         "and cases on the FDA deleted-case list were removed.")

caption("Table S1.2.", "Psychiatric event groups defined from MedDRA preferred terms (PTs).")
ev = pd.DataFrame({
    "Event group": ["Suicidality", "Depression", "Anxiety", "Eating disorders", "Psychotic symptoms"],
    "MedDRA preferred terms": [
        "Suicidal ideation; Suicide attempt; Completed suicide; Intentional self-injury; Self-injurious ideation; Suicidal behaviour",
        "Depression; Depressed mood; Major depression; Depressive symptom",
        "Anxiety; Anxiety disorder; Panic attack; Panic disorder",
        "Anorexia nervosa; Bulimia nervosa; Binge eating disorder; Eating disorder; Avoidant/restrictive food intake disorder",
        "Psychotic disorder; Hallucination; Delusion; Psychotic symptom"]})
add_table(ev)
footnote("PT matching was case-insensitive on the REAC table PT field.")

caption("Table S1.3.", "Case-level detail of the 20 pediatric GLP-1 RA suicidality reports.")
cases = pd.read_csv(f"{RES}/faers_full_ped_suicidality_cases.csv")
cases["QC flag"] = np.where(cases.age < 2, "Implausible age", "")
cases["Age"] = cases.age.map(lambda x: f"{x:.0f}" if x >= 2 else f"{x:.1f}")
cases["Seriousness"] = cases.apply(lambda r: "; ".join([s for s, c in
    [("Death", "death"), ("Life-threatening", "life_threat"), ("Hospitalization", "hospitalization")] if r[c] == 1]
    or (["Non-serious"] if r.serious == 0 else ["Serious, other"])), axis=1)
c1 = cases.rename(columns={"primaryid": "FAERS case ID", "year": "Year", "sex": "Sex",
                           "drugs": "Drug(s)", "indication_pt": "Indication PT(s)"})
c1["Indication PT(s)"] = c1["Indication PT(s)"].str.replace("|", "; ", regex=False).str.title()
c1 = c1[["FAERS case ID", "Year", "Age", "Sex", "Drug(s)", "Indication PT(s)", "Seriousness", "QC flag"]]
add_table(c1)
footnote("A sensitivity analysis excluding the implausible-age case and one suspected duplicate cluster gave "
         "a class-level pediatric suicidality ROR of 1.65 (95% CI 1.04\u20132.61). "
         "FAERS case IDs are public identifiers from the FDA FAERS quarterly data files.")

# ================= S2 =================
doc.add_page_break()
heading("S2. FAERS full disproportionality results")

def ror_table(path, strata=None):
    d = pd.read_csv(path)
    if strata: d = d[d.stratum == strata]
    out = pd.DataFrame({
        "Drug": d.drug, "Event group": d.event_group,
        "a": d.a.astype(int),
        "ROR (95% CI)": [f_or(r, l, h) for r, l, h in zip(d.ROR, d.ROR_lo, d.ROR_hi)],
        "IC025": d.IC025.map(lambda x: f"{x:.2f}"),
        "Haldane": d.haldane.map({True: "Y", False: ""}),
        "Signal": d.apply(lambda r: "ROR+IC025" if (r.signal_ROR and r.signal_IC025) else ("ROR" if r.signal_ROR else ""), axis=1)})
    return out

caption("Table S2.1.", "Pediatric (<18 years) disproportionality: ROR and Bayesian IC025 by drug and event group.")
add_table(ror_table(f"{RES}/faers_full_ror_pediatric.csv"))
footnote("a = exposed cases with the event; ROR = reporting odds ratio with 95% CI (Haldane\u2013Anscombe correction when any cell = 0, flagged Y); "
         "IC025 = lower bound of the 95% credibility interval of the WHO-UMC information component. "
         "Signal: ROR = lower 95% CI >1; ROR+IC025 = both criteria met.")

caption("Table S2.2.", "Adult (\u226518 years) disproportionality: ROR and IC025 by drug and event group.")
add_table(ror_table(f"{RES}/faers_full_ror_adult.csv"))
footnote("Conventions as in Table S4.")

per = pd.read_csv(f"{RES}/faers_full_ror_period.csv")
per = per[per.drug == "GLP1-class(all)"]
p1 = pd.DataFrame({
    "Period": per.period, "Event group": per.event_group,
    "N (pediatric GLP-1 reports)": per.N_drug.astype(int), "a": per.a.astype(int),
    "ROR (95% CI)": [f_or(r, l, h) for r, l, h in zip(per.ROR, per.ROR_lo, per.ROR_hi)],
    "IC025": per.IC025.map(lambda x: f"{x:.2f}"),
    "Signal": per.signal_ROR.map({True: "ROR", False: ""})})
caption("Table S2.3.", "Pediatric class-level disproportionality before and after 2023.")
add_table(p1)
footnote("Pre-2023: 2012Q4\u20132022Q4; post-2023: 2023Q1\u20132026Q2. Conventions as in Table S2.1.")

caption("Table S2.4.", "Indication groups of pediatric GLP-1 RA reports (n=540), overall and by drug.")
ind = pd.read_csv(f"{RES}/faers_full_ped_glp1_indication.csv")
tr = {"特殊暴露途径(自我伤害用药/宫内暴露)": "Special exposure route (self-harm use / in utero exposure)",
      "遗传/综合征性肥胖": "Genetic/syndromic obesity",
      "体重/肥胖管理": "Weight/obesity management",
      "糖尿病及糖代谢": "Diabetes and glucose metabolism",
      "其他代谢(PCOS/脂肪肝等)": "Other metabolic (PCOS, NAFLD, etc.)",
      "其他": "Other",
      "未报告/未知": "Not reported/unknown"}
ind["Indication group"] = ind.group_cn.map(tr)
ind = ind[["Indication group", "total", "Semaglutide", "Liraglutide", "Tirzepatide",
           "Exenatide", "Dulaglutide", "Lixisenatide", "Setmelanotide"]]
ind = ind.rename(columns={"total": "All"})
add_table(ind)
footnote("Indication groups were assigned from INDI table PTs; reports may carry >1 indication (assigned to a single best-fit group).")

caption("Table S2.5.", "Seriousness of pediatric GLP-1 RA reports, overall and by event group.")
ser = pd.read_csv(f"{RES}/faers_full_ped_glp1_serious.csv")
tr2 = {"GLP-1 儿科全部": "All pediatric GLP-1 RA reports",
       "其中: 自杀相关": "of which: suicidality",
       "其中: 抑郁": "of which: depression",
       "其中: 焦虑": "of which: anxiety",
       "其中: 进食障碍": "of which: eating disorders",
       "其中: 精神病性症状": "of which: psychotic symptoms"}
ser["Subset"] = ser.subset.map(tr2).fillna(ser.subset)
s1 = pd.DataFrame({
    "Subset": ser.Subset, "N": ser.N,
    "Serious, n (%)": ser.apply(lambda r: f"{int(r.serious_n)} ({r.serious_pct})", axis=1),
    "Death": ser.death, "Life-threatening": ser.life_threat, "Hospitalization": ser.hospitalization,
    "Disability": ser.disability, "Required intervention": ser.req_intervention,
    "Other serious": ser.other_serious_OT})
add_table(s1)
footnote("Seriousness categories follow FDA outcome codes and are not mutually exclusive.")

# ================= S3 =================
doc.add_page_break()
heading("S3. YRBS 2025 full model tables (survey-weighted)")
para("All models are survey-weighted logistic regressions (WT, STRATUM, PSU) adjusted for sex, "
     "race/ethnicity, and age, on the 2025 national Youth Risk Behavior Survey analytic sample.",
     size=9, space_after=8)

OC = {"持续悲伤/绝望": "Persistent sadness/hopelessness",
      "认真考虑过自杀": "Seriously considered suicide",
      "制定自杀计划": "Made a suicide plan",
      "自杀未遂": "Suicide attempt"}

prev = pd.read_csv(f"{RES}/yrbs_2025_prevalence.csv")
prev["Outcome"] = prev.outcome.map(OC)
p1 = pd.DataFrame({"BMI category": prev.bmi, "Outcome": prev.Outcome, "n": prev.n,
                   "Weighted % (95% CI)": prev.apply(lambda r: f"{r.pct:.1f} ({r.lo:.1f}\u2013{r.hi:.1f})", axis=1)})
caption("Table S3.1.", "Weighted prevalence of mental-health outcomes by BMI category.")
add_table(p1)

orr = pd.read_csv(f"{RES}/yrbs_2025_or.csv")
orr["Outcome"] = orr.outcome.map(OC)
o1 = pd.DataFrame({"Outcome": orr.Outcome, "Exposure (ref = normal weight)": orr.term,
                   "Adjusted OR (95% CI)": [f_or(r, l, h) for r, l, h in zip(orr.OR, orr.lo, orr.hi)],
                   "P": orr.p.map(f_p)})
caption("Table S3.2.", "BMI category and mental-health outcomes (reference: normal weight).")
add_table(o1)

wc = pd.read_csv(f"{RES}/yrbs_2025_weightcontrol_or.csv")
def tr_lab(s):
    for k, v in OC.items(): s = s.replace(k, v)
    s = s.replace("减重意向→", "Trying to lose weight \u2192 ").replace("感知超重(实际正常)→", "Perceived overweight despite normal weight \u2192 ")
    s = s.replace("(全体, 校正)", "(all, adjusted)").replace("限Normal", "restricted to normal weight")
    s = s.replace("限Overweight", "restricted to overweight").replace("限Obese", "restricted to obesity")
    return s
w1 = pd.DataFrame({"Analysis": wc.analysis.map(tr_lab),
                   "Adjusted OR (95% CI)": [f_or(r, l, h) for r, l, h in zip(wc.OR, wc.lo, wc.hi)],
                   "P": wc.p.map(f_p)})
caption("Table S3.3.", "Weight-loss intention and weight-perception mismatch: adjusted associations.")
add_table(w1)
footnote("Models additionally adjusted for sex, race/ethnicity, and age; stratum-restricted models drop the corresponding BMI term.")

bs = pd.read_csv(f"{RES}/yrbs_2025_bysex.csv")
bs["Outcome"] = bs.outcome.map(OC)
b1 = pd.DataFrame({"Sex": bs.sex, "BMI category": bs.bmi, "Outcome": bs.Outcome, "n": bs.n,
                   "Weighted % (95% CI)": bs.apply(lambda r: f"{r.pct:.1f} ({r.lo:.1f}\u2013{r.hi:.1f})", axis=1)})
caption("Table S3.4a.", "Outcome prevalence by sex and BMI category.")
add_table(b1)

tl = pd.read_csv(f"{RES}/yrbs_2025_trylose_bysex.csv")
t1 = pd.DataFrame({"Sex": tl.sex, "BMI category": tl.bmi, "n": tl.n,
                   "Trying to lose weight, weighted % (95% CI)": tl.apply(lambda r: f"{min(r.pct,99.9):.1f} ({max(r.lo,0):.1f}\u2013{min(r.hi,100):.1f})", axis=1)})
caption("Table S3.4b.", "Prevalence of trying to lose weight by sex and BMI category.")
add_table(t1)
footnote("CIs truncated at 0\u2013100% for presentation.")

# ================= S4 =================
doc.add_page_break()
heading("S4. Mendelian randomization supporting analyses")
para("Instrument selection: variants within \u00b1100 kb of GLP1R (or comparator loci), P<5\u00d710\u207b\u2078, "
     "clumped with PLINK (r\u00b2=0.3, 100 kb, 1000 Genomes EUR). Estimates are fixed-effect IVW "
     "(Wald ratio when n=1) per 1-SD higher BMI z-score; the BMI-lowering (drug-mimicking) direction is the reciprocal.",
     size=9, space_after=8)

# S4.1 instruments
inst_frames = []
wave_sets = {"3months": ["rs11963172", "rs116623469", "rs9470976", "rs2268657"],
             "6months": ["rs11963172", "rs1820722"],
             "8months": ["rs11963172", "rs1820721"],
             "1year": ["rs1412265", "rs28360623", "rs179269"],
             "1.5years": ["rs1412265", "rs1820721"]}
wl = {"3months": "3 months", "6months": "6 months", "8months": "8 months",
      "1year": "1 year", "1.5years": "1.5 years"}
for w, snps in wave_sets.items():
    d = pd.read_csv(f"{MR}/instruments_MoBa2022_{w}_GLP1R_100kb.csv")
    d["rs_key"] = d.rsid.str.split("|").str[0]
    d = d[d.rs_key.isin(snps)]
    inst_frames.append(pd.DataFrame({
        "Wave": wl[w], "rsID": d.rsid, "Position (GRCh37)": d.pos.map(lambda x: f"{int(x):,}"),
        "EA/NEA": d.ea + "/" + d.nea, "EAF": d.eaf.map(lambda x: f"{x:.3f}"),
        "\u03b2 (SE)": d.apply(lambda r: f"{r.beta:.4f} ({r.se:.4f})", axis=1),
        "P": d.p.map(lambda x: f"{x:.1e}"), "F": d.F.map(lambda x: f"{x:.1f}")}))
i1 = pd.concat(inst_frames)
caption("Table S4.1.", "Wave-specific GLP1R cis-instruments for childhood BMI z-score (MoBa).")
add_table(i1)
footnote("EA = effect (BMI-increasing) allele; NEA = other allele; EAF = effect-allele frequency; "
         "rs2268657 was additionally matched to the exome-chip probe exm2257793 in outcome statistics (position-based fallback).")

# S4.2 wave-specific results
res = pd.read_csv(f"{MR}/child_mr_results.csv")
g = res[(res.gene == "GLP1R") & res.dataset.str.startswith("MoBa2022")].copy()
g["Wave"] = g.dataset.str.replace("MoBa2022_", "").map(wl)
g = g.sort_values(["Wave", "model"])
r1 = pd.DataFrame({"Wave": g.Wave, "Model": g.model.str.replace("Model", "Model "),
                   "n SNP": g.n_snp,
                   "OR per 1-SD higher BMI (95% CI)": [f_or(r, l, h) for r, l, h in zip(g.or_, g.or_lo, g.or_hi)],
                   "P": g.p.map(f_p), "P (heterogeneity)": g.p_het.map(f_p),
                   "Steiger pass": g.steiger_all_ok.map({True: "Y", False: "N"})})
caption("Table S4.2.", "Wave-specific GLP1R cis-MR estimates, iPSYCH Models 1\u20133 (suicide attempt).")
add_table(r1)
footnote("Model 1 = any suicide attempt; Model 2 = attempt requiring hospital contact; Model 3 = attempt within affective-disorder subgroup.")

# S4.3 developmental scan
dev = pd.read_csv(f"{MR}/child_mr_developmental.csv")
dev = dev[dev.gene == "GLP1R"].copy()
wave_order = {w: i for i, w in enumerate(dev.timepoint.unique())}
dev["o"] = dev.timepoint.map(wave_order)
dev = dev.sort_values(["o", "model"])
d1 = pd.DataFrame({"Wave": dev.timepoint, "Age (y)": dev.age, "Model": dev.model.str.replace("Model", "M"),
                   "n SNP": dev.n_snp, "Mean F": dev.mean_F.map(lambda x: f"{x:.1f}"),
                   "OR per 1-SD higher BMI (95% CI)": [f_or(r, l, h) for r, l, h in zip(dev.or_, dev.or_lo, dev.or_hi)],
                   "P": dev.p.map(f_p),
                   "OR per 1-SD lower BMI (95% CI)": [f_or(r, l, h) for r, l, h in zip(dev.or_down, dev.or_down_lo, dev.or_down_hi)]})
caption("Table S4.3.", "Fixed infant-instrument scan across twelve MoBa waves (birth to 8 years).")
add_table(d1)
footnote("The 9-variant fixed infant set was carried forward across waves; waves with mean F<10 (5\u20138 years) are shown for illustration only.")

# S4.4 comparators
c = res[res.gene.isin(["MC4R", "GIPR"])].copy()
c1 = pd.DataFrame({"Locus": c.gene, "Exposure": c.dataset, "Window": c.window,
                   "Model": c.model.str.replace("Model", "M"), "n SNP": c.n_snp,
                   "OR per 1-SD higher BMI (95% CI)": [f_or(r, l, h) for r, l, h in zip(c.or_, c.or_lo, c.or_hi)],
                   "P": c.p.map(f_p)})
caption("Table S4.4.", "Comparator loci (MC4R, GIPR) cis-MR estimates.")
add_table(c1)

# S4.5 adult
ad = pd.read_csv(f"{MR}/adult_meta_mvp_mr.csv")
a1 = pd.DataFrame({"Exposure": ad.dataset, "Model": ad.model.str.replace("Model", "M"),
                   "n harmonized": ad.n_harm, "Mean F": ad.mean_F.map(lambda x: f"{x:.1f}"),
                   "OR per 1-SD higher BMI (95% CI)": [f_or(r, l, h) for r, l, h in zip(ad.or_, ad.or_lo, ad.or_hi)],
                   "P": ad.p.map(f_p), "P (het)": ad.p_het.map(f_p)})
caption("Table S4.5a.", "Adult GLP1R cis-MR: GIANT\u2013Million Veteran Program meta-analysis (N\u22481.11M) and MVP-only.")
add_table(a1)

rx = pd.read_csv(f"{MR}/adult_glp1r_relaxed.csv")
rx1 = rx.dropna(subset=["or_"]).copy()
r2 = pd.DataFrame({"P threshold": rx1.threshold, "Model": rx1.model.str.replace("Model", "M"),
                   "n harmonized": rx1.n_harm,
                   "OR per 1-SD higher BMI (95% CI)": [f_or(r, l, h, 2) for r, l, h in zip(rx1.or_, rx1.or_lo, rx1.or_hi)],
                   "P": rx1.p.map(f_p)})
caption("Table S4.5b.", "Relaxed-threshold sensitivity for adult BMI (ieu-b-40) instruments.")
add_table(r2)
footnote("Single-variant Wald ratios dominate at relaxed thresholds; wide CIs reflect limited instrument count.")

# S4.6 TwoSampleMR replication
ua = pd.read_csv("E:/CM/GLP1/results/moba2022_bmi_3months_Model1_20260927223339/05.Model1_mr_analysis.csv")
uo = pd.read_csv("E:/CM/GLP1/results/moba2022_bmi_3months_Model1_20260927223339/06.Model1_or_result.csv")
u = ua.merge(uo[["method", "or", "or_lci95", "or_uci95"]], on="method")
t6 = pd.DataFrame({"Method": u.method, "n SNP": u.nsnp,
                   "\u03b2 (SE)": u.apply(lambda r: f"{r.b:.3f} ({r.se:.3f})", axis=1),
                   "OR (95% CI)": [f_or(r, l, h) for r, l, h in zip(u["or"], u.or_lci95, u.or_uci95)],
                   "P": u.pval.map(f_p)})
caption("Table S4.6.", "Independent replication with TwoSampleMR (R) at the 3-month wave (3-SNP set).")
add_table(t6)
footnote("Heterogeneity: IVW Q=0.55, df=2, P=.759; MR-Egger Q=0.43, df=1, P=.514. "
         "Our Python IVW on the 4-variant wave set (adding rs2268657 via exome-probe fallback) gave OR 0.55 (0.39\u20130.79), P=.0013.")

# S4.7 GLS
gl = pd.read_csv(f"{MR}/ld_gls_sensitivity.csv")
g7 = pd.DataFrame({"Wave": gl.wave.map(wl), "n SNP": gl.n_snp,
                   "GLS OR per 1-SD higher BMI (95% CI)": [f_or(r, l, h) for r, l, h in zip(gl.or_gls, gl.or_lo, gl.or_hi)],
                   "P": gl.p_gls.map(f_p),
                   "GLS OR per 1-SD lower BMI (95% CI)": [f_or(r, l, h) for r, l, h in zip(gl.or_down, gl.or_down_lo, gl.or_down_hi)]})
caption("Table S4.7.", "LD-corrected generalized-least-squares sensitivity (Model 1).")
add_table(g7)
footnote("GLS uses the 1000 Genomes EUR LD correlation matrix among instruments (PLINK --r square); "
         "exm2257793 is absent from the panel and the 3-month GLS uses the remaining 4 variants. "
         "Estimates are essentially unchanged from independence-assumption IVW.")

# S4.8 LD matrix
M = pd.read_csv(f"{MR}/ld_matrix_infant_adult_r2.csv", index_col=0).round(2)
M.insert(0, "Infant variant", M.index)
caption("Table S4.8.", "Pairwise LD (r\u00b2, 1000 Genomes EUR) between infant-wave instruments (rows) and adult harmonized variants (columns).")
add_table(M, font=7.5)
footnote("Adult variants: meta-analysis harmonized set (rs12213929, rs2268647, rs9470954) plus MVP-only set. "
         "Note rs1412265 (1\u20131.5-year waves) is in strong LD with the adult lead rs12213929 (r\u00b2=0.97); "
         "the 3-month lead cluster (rs9470976/rs11963172) is largely independent of adult leads (r\u00b2\u22640.12).")

# ================= S5 =================
doc.add_page_break()
heading("S5. Figure S1")
p = doc.add_paragraph()
r = p.add_run("Figure S1. "); r.bold = True; r.font.size = Pt(9)
r2 = p.add_run("Wave-specific GLP1R cis-MR estimates (Models 1\u20133) for suicide attempt per 1-SD higher "
               "childhood BMI z-score. Fixed-effect IVW (Wald ratio for single-variant sets); error bars, 95% CI. "
               "Moved from the main text for compactness.")
r2.font.size = Pt(9)
doc.add_picture(f"{MR}/Fig2_forest_pub.png", width=Inches(6.7))

doc.save("manuscript/Supplementary_Appendix_v2.docx")
print("saved manuscript/Supplementary_Appendix_v2.docx")
