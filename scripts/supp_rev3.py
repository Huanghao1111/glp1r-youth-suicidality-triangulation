#!/usr/bin/env python3
"""Round-2 supplement edits:
1. Table S1: fix adult/age-missing N (canonical dedup view: 9,457,686 / 7,511,216).
2. Table S25: rebuild as six-class comparator panel + GI control rows (15 rows).
3. Update S25 caption (para) and footnote.
4. Insert Table S28 (coloc posteriors) after S27 footnote, before Figure S1 caption.
"""
import copy
import docx

path = 'manuscript/submission_package/Supplementary_Appendix.docx'
s = docx.Document(path)

def set_cell(cell, text):
    """replace text preserving the first run's formatting"""
    p = cell.paragraphs[0]
    runs = p.runs
    if not runs:
        p.add_run(text)
        return
    runs[0].text = text
    for r in runs[1:]:
        r.text = ''

# ---------- 1. Table S1 N fix ----------
t1 = s.tables[0]
nfix = 0
for row in t1.rows:
    for cell in row.cells:
        for p in cell.paragraphs:
            for r in p.runs:
                if '9,457,685' in r.text:
                    r.text = r.text.replace('9,457,685', '9,457,686'); nfix += 1
                if '7,511,217' in r.text:
                    r.text = r.text.replace('7,511,217', '7,511,216'); nfix += 1
print('S1 cells fixed:', nfix)

# ---------- 2. rebuild Table S25 (docx table 24) ----------
t = s.tables[24]
ROWS = [
    ["Drug class", "Period", "Exposed, n", "Suicidality cases, a", "ROR (95% CI)", "IC025", "Fold (post/pre)"],
    ["GLP-1 RA (index)", "pre-2023", "130", "3", "1.01 (0.32–3.17)", "−2.23", ""],
    ["GLP-1 RA (index)", "post-2023", "410", "17", "2.21 (1.36–3.60)", "0.21", "2.19"],
    ["SSRIs (stable-warning comparator)", "pre-2023", "10,538", "1,959", "12.06 (11.42–12.73)", "2.95", ""],
    ["SSRIs (stable-warning comparator)", "post-2023", "3,824", "841", "17.94 (16.49–19.52)", "3.40", "1.49"],
    ["Second-generation antipsychotics (stable-warning comparator)", "pre-2023", "17,116", "1,494", "4.67 (4.41–4.94)", "1.85", ""],
    ["Second-generation antipsychotics (stable-warning comparator)", "post-2023", "6,027", "608", "6.57 (6.00–7.19)", "2.26", "1.41"],
    ["GLP-1 RA — GI control events", "pre-2023", "130", "32", "4.59 (3.08–6.84)", "1.21", ""],
    ["GLP-1 RA — GI control events", "post-2023", "410", "106", "5.18 (4.15–6.47)", "1.69", "1.13"],
]
# existing table: 11 rows x 7 cols -> need 10 rows; rebuild then drop trailing rows
while len(t.rows) < len(ROWS):
    t.add_row()
for ri, vals in enumerate(ROWS):
    for ci, v in enumerate(vals):
        set_cell(t.rows[ri].cells[ci], v)
while len(t.rows) > len(ROWS):
    t._tbl.remove(t.rows[-1]._tr)
print('S25 rebuilt:', len(t.rows), 'rows')

# ---------- 3. caption + footnote ----------
for r in s.paragraphs[54].runs:
    if "ADHD comparator class and gastrointestinal control events" in r.text:
        r.text = r.text.replace("ADHD comparator class and gastrointestinal control events",
                                "stable-warning psychotropic comparator classes and gastrointestinal control events")
print('S25 caption:', s.paragraphs[54].text[:120])

# footnote para 55: full replacement (keep formatting of first run)
p55 = s.paragraphs[55]
new_note = ("Comparator classes were chosen a priori as psychotropic drug classes with established, long-standing "
            "suicidality warnings and no 2023 class-level media event, providing the reference background for pediatric "
            "suicidality reporting intensity (primary/secondary suspect, mirroring the GLP-1 RA exposure definition): "
            "SSRIs (fluoxetine, sertraline, escitalopram, citalopram, fluvoxamine, paroxetine) and second-generation "
            "antipsychotics (aripiprazole, risperidone, quetiapine, olanzapine, ziprasidone, lurasidone, paliperidone). "
            "GI control events: nausea, vomiting, diarrhoea, constipation, pancreatitis — expected GLP-1 RA pharmacology "
            "without suicidality-related media attention. Against the stable-warning classes' 1.4- to 1.5-fold post-2023 "
            "rise, the GLP-1 RA rise (2.19-fold) leaves a modest excess; whether it is drug-specific or reflects the "
            "class's intense media attention cannot be separated in this source. Fold = post-2023 ROR / pre-2023 ROR.")
if p55.runs:
    p55.runs[0].text = new_note
    for r in p55.runs[1:]:
        r.text = ''
print('S25 footnote updated')

# ---------- 4. insert Table S28 (coloc) after para 59 (S27 footnote) ----------
COLOC = [
    ["Exposure", "Outcome", "Variants, n", "PP.H3 (distinct)", "PP.H4 (shared)", "Top shared variant", "Top P (exposure / outcome)"],
    ["MoBa 3-month BMI z", "iPSYCH Model 1", "1,901", "0.076", "0.165", "rs9470972", "7.7×10⁻¹⁰ / .012"],
    ["MoBa 3-month BMI z", "iPSYCH Model 2", "1,901", "0.085", "0.092", "rs9470972", "7.7×10⁻¹⁰ / .026"],
    ["MoBa 3-month BMI z", "iPSYCH Model 3", "1,900", "0.107", "0.083", "rs11963172", "2.2×10⁻⁹ / .028"],
    ["MoBa 1-year BMI z", "iPSYCH Model 1", "1,901", "0.082", "0.096", "rs179269", "2.7×10⁻⁸ / .008"],
    ["Adult BMI meta-analysis", "iPSYCH Model 1", "2,124", "0.100", "0.018", "rs12213929", "2.0×10⁻¹⁵ / .23"],
    ["Adult BMI meta-analysis", "iPSYCH Model 3", "2,124", "0.128", "0.034", "rs12213929", "2.0×10⁻¹⁵ / .13"],
]
cap_template = s.paragraphs[58]  # "Table S27. ..." caption
note_template = s.paragraphs[59] # S27 footnote
anchor = s.paragraphs[59]._p     # insert after this

# caption paragraph
cap = copy.deepcopy(cap_template._p)
anchor.addnext(cap)
from docx.text.paragraph import Paragraph
cap_p = Paragraph(cap, s.paragraphs[58]._parent)
for r in cap_p.runs[1:]:
    r.text = ''
cap_p.runs[0].text = "Table S28. GLP1R colocalization posteriors (coloc.abf), GLP1R ±250-kb region."

# table: deepcopy S25 (7 cols, 15 rows) as structure template, shrink to 7 rows
tbl_xml = copy.deepcopy(s.tables[24]._tbl)
cap.addnext(tbl_xml)
from docx.table import Table
t28 = Table(tbl_xml, s.tables[24]._parent)
while len(t28.rows) > len(COLOC):
    t28._tbl.remove(t28.rows[-1]._tr)
for ri, vals in enumerate(COLOC):
    for ci, v in enumerate(vals):
        set_cell(t28.rows[ri].cells[ci], v)

# footnote paragraph after the table
note = copy.deepcopy(note_template._p)
tbl_xml.addnext(note)
note_p = Paragraph(note, s.paragraphs[59]._parent)
for r in note_p.runs[1:]:
    r.text = ''
note_p.runs[0].text = ("coloc.abf (Wakefield approximate Bayes factors; priors p1=p2=10⁻⁴, p12=10⁻⁵; W=0.15 quantitative, "
    "0.20 case–control), GRCh37. PP.H0–H2 ≈ 0/0.76–0.88/≈0 (exposure-only hypothesis H1 dominates) and are omitted for "
    "brevity. No suicide-attempt variant in the region approaches genome-wide significance (top outcome P=.008), so the "
    "analysis is underpowered to distinguish distinct (H3) from shared (H4) causal variants; posteriors are uninformative "
    "rather than evidence against colocalization. The single-causal-variant assumption is additionally strained by 2–3 "
    "near-independent BMI signals in the region.")

# renumber any duplicate paraIds introduced by deepcopy
import re as _re
seen = set()
for p in s.paragraphs:
    pid = p._p.get('{http://schemas.openxmlformats.org/wordprocessingml/2014/main}paraId')
    if pid:
        if pid in seen:
            import random
            p._p.set('{http://schemas.openxmlformats.org/wordprocessingml/2014/main}paraId',
                     f'{random.getrandbits(32):08X}')
        else:
            seen.add(pid)

s.save(path)
print('S28 inserted; saved')

# verify
s2 = docx.Document(path)
print('tables now:', len(s2.tables))
for i,p in enumerate(s2.paragraphs):
    if 'Table S28' in p.text or 'coloc.abf' in p.text:
        print(f'[{i}] {p.text[:90]}')
EOF
