# -*- coding: utf-8 -*-
"""LD-corrected GLS sensitivity for wave-specific GLP1R instruments (Model 1)
+ pairwise LD matrix between infant and adult lead variants.
PLINK 1.9 EUR panel at E:/CM/GLP1/LD/bfile/EUR."""
import os, subprocess
import numpy as np
import pandas as pd

WS = "results/mr_cis"
DIR = "E:/CM/GLP1/mr_input"
PLINK = "E:/CM/GLP1/LD/plink_win64/plink.exe"
BFILE = "E:/CM/GLP1/LD/bfile/EUR"

# --- collect wave-specific instrument sets (Model 1, GLP1R MoBa waves) ---
res = pd.read_csv(f"{WS}/child_mr_results.csv")
res = res[(res.gene == "GLP1R") & res.dataset.str.startswith("MoBa2022") & (res.model == "Model1")]
waves = {}
for _, r in res.iterrows():
    snps = []
    for s in r.snps.split(","):
        snps.extend(s.split("|"))  # exome probe fallbacks
    waves[r.dataset.replace("MoBa2022_", "")] = snps
print({k: v for k, v in waves.items()})

# adult harmonized SNPs
meta = pd.read_csv(f"{DIR}/mrinput_AdultBMI_meta_GLP1R_100kb_vs_iPSYCH_Model1.csv")
mvp = pd.read_csv(f"{DIR}/mrinput_AdultBMI_MVP_GLP1R_100kb_vs_iPSYCH_Model1.csv")
adult_meta = list(meta.rsid); adult_mvp = list(mvp.rsid)
print("adult meta:", adult_meta); print("adult mvp:", adult_mvp)

infant_union = sorted({s for v in waves.values() for s in v})
all_snps = sorted(set(infant_union) | set(adult_meta) | set(adult_mvp))
print("union n =", len(all_snps))

# --- PLINK pairwise r2 (square, signed r via --r yes-really?) ---
lst = f"{DIR}/ld_snplist.txt"
with open(lst, "w") as f:
    f.write("\n".join(all_snps))
outpref = f"{DIR}/ld_glp1r_all"
subprocess.run([PLINK, "--bfile", BFILE, "--extract", lst, "--r", "square",
                "--out", outpref, "--silent"], capture_output=True, text=True)
# plink --r square writes .ld (square matrix) and .nosex; variant order = .ld rows
ldf = f"{outpref}.ld"
R = pd.read_csv(ldf, sep=r"\s+", header=None, engine="python")
# panel variants actually present
import io
found = subprocess.run([PLINK, "--bfile", BFILE, "--extract", lst, "--write-snplist",
                        "--out", f"{DIR}/ld_present", "--silent"], capture_output=True, text=True)
present = [l.strip() for l in open(f"{DIR}/ld_present.snplist") if l.strip()]
print("present in EUR panel:", len(present), "of", len(all_snps))
missing = [s for s in all_snps if s not in present]
print("missing:", missing)
R.index = present; R.columns = present

# --- GLS per wave (Model 1), fixed effect ---
rows = []
for wave, snps in waves.items():
    mi = pd.read_csv(f"{DIR}/mrinput_MoBa2022_{wave}_GLP1R_100kb_vs_iPSYCH_Model1.csv")
    mi["rs_key"] = mi.rsid.str.split("|").str[0]
    mi = mi[mi.rs_key.isin(present)]
    if len(mi) < 2:
        # try fallback IDs
        mi = pd.read_csv(f"{DIR}/mrinput_MoBa2022_{wave}_GLP1R_100kb_vs_iPSYCH_Model1.csv")
        mi["rs_key"] = mi.rsid.str.split("|").str[-1]
        mi = mi[mi.rs_key.isin(present)]
    if len(mi) < 2:
        print(wave, "too few SNPs with LD:", len(mi)); continue
    x = mi.beta.values.astype(float)
    y = mi.beta_out_aligned.values.astype(float)
    sey = mi.se_out.values.astype(float)
    Rw = R.loc[mi.rs_key, mi.rs_key].values.astype(float)
    Om = Rw * np.outer(sey, sey)
    Om += np.eye(len(mi)) * 1e-12
    Oi = np.linalg.inv(Om)
    beta = (x @ Oi @ y) / (x @ Oi @ x)
    se = np.sqrt(1.0 / (x @ Oi @ x))
    z = beta / se
    from math import erfc
    p = erfc(abs(z) / np.sqrt(2))
    rows.append(dict(wave=wave, n_snp=len(mi), dropped=",".join(
        sorted(set(snps) - set(mi.rs_key))), beta_gls=beta, se_gls=se,
        or_gls=np.exp(beta), or_lo=np.exp(beta - 1.96 * se), or_hi=np.exp(beta + 1.96 * se),
        p_gls=p, or_down=1/np.exp(beta), or_down_lo=1/np.exp(beta + 1.96 * se),
        or_down_hi=1/np.exp(beta - 1.96 * se)))
    print(wave, f"n={len(mi)} GLS OR={np.exp(beta):.3f} ({np.exp(beta-1.96*se):.3f}-{np.exp(beta+1.96*se):.3f}) p={p:.2e}")

gls = pd.DataFrame(rows)
gls.to_csv(f"{WS}/ld_gls_sensitivity.csv", index=False)

# --- infant x adult LD matrix (r2) ---
sub_i = [s for s in infant_union if s in present]
sub_a = sorted(set(adult_meta) | set(adult_mvp))
sub_a = [s for s in sub_a if s in present]
M = (R.loc[sub_i, sub_a] ** 2)
M.to_csv(f"{WS}/ld_matrix_infant_adult_r2.csv")
print("LD matrix saved:", M.shape, "max r2 =", M.values.max())
