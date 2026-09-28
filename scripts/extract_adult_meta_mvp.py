#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Extract GLP1R cis windows from user's BMI_meta.txt (ieu-b-40 + MVP meta, N~1.1M)
and BMI_MVP_European.txt (MVP European, N=424k).
MVP/meta are GRCh38; shift vs hg19 derived empirically from shared rsIDs, then all
positions are stored hg19-mapped for consistency with the rest of the pipeline.
Outputs: AdultBMI_meta_GLP1R_window.csv, AdultBMI_MVP_GLP1R_window.csv
(columns match extract_cis.py convention)."""
import pandas as pd
import numpy as np

DIR = "E:/CM/GLP1/mr_input"
G19_START, G19_END = 39_016_574, 39_055_519   # GLP1R hg19
W250, W100 = 250_000, 100_000

# --- hg19 reference positions from existing ieu-b-40 window table ---
ref = pd.read_csv(f"{DIR}/AdultBMI_ieub40_GLP1R_window.csv")
ref = ref[["rsid", "pos", "in_100kb"]].drop_duplicates("rsid")
pos19 = dict(zip(ref["rsid"], ref["pos"]))
print(f"hg19 reference rsIDs in +-250kb window: {len(ref)}")

# --- pass 1: scan MVP (has chr/pos, GRCh38) ---
keep_rows = []
for chunk in pd.read_csv("E:/CM/GLP1/BMI_MVP_European.txt", sep="\t", chunksize=2_000_000,
                         dtype={"chr": str}):
    sub = chunk[(chunk["chr"] == "6") & (chunk["pos"] >= 38_500_000) & (chunk["pos"] <= 39_700_000)]
    if len(sub):
        keep_rows.append(sub)
mvp = pd.concat(keep_rows, ignore_index=True)
print(f"MVP chr6 38.5-39.7Mb(GRCh38) rows: {len(mvp)}")

# --- empirical GRCh38->hg19 shift from shared rsIDs ---
both = mvp[mvp["SNP"].isin(pos19)].copy()
both["pos19"] = both["SNP"].map(pos19)
shifts = both["pos"] - both["pos19"]
SHIFT = int(shifts.mode().iloc[0])
agree = (shifts == SHIFT).mean()
print(f"shift GRCh38-hg19 = {SHIFT:+d} bp, agreement {agree:.1%} of {len(shifts)} shared rsIDs")
assert agree > 0.95, "shift not constant - check build!"

# --- window membership on hg19-mapped coordinates ---
mvp["pos19"] = mvp["pos"] - SHIFT
m_in = mvp[(mvp["pos19"] >= G19_START - W250) & (mvp["pos19"] <= G19_END + W250)].copy()
m_in["in_100kb"] = (m_in["pos19"] >= G19_START - W100) & (m_in["pos19"] <= G19_END + W100)
m_in["dist_to_gene"] = np.where(m_in["pos19"] < G19_START, m_in["pos19"] - G19_START,
                                np.where(m_in["pos19"] > G19_END, m_in["pos19"] - G19_END, 0))
m_out = pd.DataFrame(dict(
    rsid=m_in["SNP"], chr=6, pos=m_in["pos19"], pos38=m_in["pos"],
    ea=m_in["effect_allele"], nea=m_in["other_allele"], eaf=m_in["eaf"],
    info=np.nan, beta=m_in["beta"], se=m_in["se"], p=m_in["pval"], n=m_in["samplesize"],
    dataset="AdultBMI_MVP", gene="GLP1R", dist_to_gene=m_in["dist_to_gene"],
    in_100kb=m_in["in_100kb"]))
m_out["F"] = (m_out["beta"] / m_out["se"]) ** 2
m_out.to_csv(f"{DIR}/AdultBMI_MVP_GLP1R_window.csv", index=False)
c5 = m_out[m_out["in_100kb"] & (m_out["p"] < 5e-8)]
print(f"MVP window: {len(m_out)} SNPs (+-250kb), {m_out['in_100kb'].sum()} in +-100kb, "
      f"{len(c5)} pass p<5e-8, min p={m_out[m_out['in_100kb']]['p'].min():.2e}")

# --- pass 2: scan meta (rsID only, no positions) ---
rsid_set = set(pos19) | set(m_out["rsid"])
meta_rows = []
for chunk in pd.read_csv("E:/CM/GLP1/BMI_meta.txt", sep="\t", chunksize=2_000_000):
    sub = chunk[chunk["SNP"].isin(rsid_set)]
    if len(sub):
        meta_rows.append(sub)
meta = pd.concat(meta_rows, ignore_index=True)
print(f"\nmeta rows matching window rsID set: {len(meta)}")
meta["pos19"] = meta["SNP"].map(pos19)
# MVP-only rsIDs: fill pos19 from the MVP extraction
mvp_pos = dict(zip(m_out["rsid"], m_out["pos"]))
miss = meta["pos19"].isna()
meta.loc[miss, "pos19"] = meta.loc[miss, "SNP"].map(mvp_pos)
meta = meta.dropna(subset=["pos19"])
meta["pos19"] = meta["pos19"].astype(int)
meta["in_100kb"] = (meta["pos19"] >= G19_START - W100) & (meta["pos19"] <= G19_END + W100)
meta["dist_to_gene"] = np.where(meta["pos19"] < G19_START, meta["pos19"] - G19_START,
                                np.where(meta["pos19"] > G19_END, meta["pos19"] - G19_END, 0))
meta_out = pd.DataFrame(dict(
    rsid=meta["SNP"], chr=6, pos=meta["pos19"],
    ea=meta["effect_allele"], nea=meta["other_allele"], eaf=meta["eaf"],
    info=meta["HetISq"], beta=meta["beta"], se=meta["se"], p=meta["pval"], n=meta["samplesize"],
    dataset="AdultBMI_meta", gene="GLP1R", dist_to_gene=meta["dist_to_gene"],
    in_100kb=meta["in_100kb"]))
meta_out["F"] = (meta_out["beta"] / meta_out["se"]) ** 2
meta_out.to_csv(f"{DIR}/AdultBMI_meta_GLP1R_window.csv", index=False)
c5m = meta_out[meta_out["in_100kb"] & (meta_out["p"] < 5e-8)]
print(f"meta window: {len(meta_out)} SNPs, {meta_out['in_100kb'].sum()} in +-100kb, "
      f"{len(c5m)} pass p<5e-8, min p={meta_out[meta_out['in_100kb']]['p'].min():.2e}")
print("\ntop meta SNPs (+-100kb):")
top = meta_out[meta_out["in_100kb"]].nsmallest(12, "p")
print(top[["rsid", "pos", "ea", "nea", "beta", "se", "p", "n", "F", "info"]].to_string(index=False))
