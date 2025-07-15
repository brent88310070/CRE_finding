#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
merge_tomato_deg.py
===================
Merge per-sample DEG result TSVs, then output stringent DEG / Non-DEG lists.

每個輸入 TSV 必含欄位：Geneid, log2FoldChange, lfcSE, padj
"""

from pathlib import Path
from functools import partial
import argparse
import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues


# ----------------------------------------------------------------------
# Argument parser
# ----------------------------------------------------------------------
def get_args():
    p = argparse.ArgumentParser(
        description="Merge multiple DEG tables and call rigorous "
                    "DEG / Non-DEG gene sets"
    )
    p.add_argument("--input_dir", default="./tomato_deg_results",
                   help="Folder containing per-sample *.tsv DEG tables")
    p.add_argument("--out_dir",   default="./deg_summary",
                   help="Folder to write merged summary files")
    p.add_argument("--padj_th",   type=float, default=0.05,
                   help="FDR threshold for significance")
    p.add_argument("--fc_th",     type=float, default=1.0,
                   help="abs(log2FC) threshold for significance")
    p.add_argument("--prefix",    default="tomato",
                   help="Prefix for output filenames (e.g. species)")
    return p.parse_args()


# ----------------------------------------------------------------------
def fisher_p(row, sig_cols, min_p=1e-300):
    """Fisher’s method to combine p-values along a row."""
    ps = row[sig_cols].dropna()
    if ps.empty:
        return np.nan
    ps = np.where(ps == 0, min_p, ps)
    return combine_pvalues(ps)[1]


def meta_fc(row, fc_cols, se_cols):
    """Inverse-variance weighted mean log2FC."""
    fcs = row[fc_cols].values
    ses = row[se_cols].values
    mask = (~np.isnan(fcs)) & (~np.isnan(ses)) & (ses > 0)
    if mask.sum() == 0:
        return np.nan
    w = 1 / (ses[mask] ** 2)
    return np.average(fcs[mask], weights=w)


# ----------------------------------------------------------------------
def main():
    args = get_args()

    INPUT_DIR  = Path(args.input_dir)
    TSV_DIR    = Path(args.out_dir)
    PADJ_TH    = args.padj_th
    FC_TH      = args.fc_th
    OUT_DEG    = TSV_DIR / f"{args.prefix}_DEG.tsv"
    OUT_NON    = TSV_DIR / f"{args.prefix}_nonDEG.tsv"

    # ---------- Merge all sample tables ----------
    dfs = []
    for path in INPUT_DIR.glob("*.tsv"):
        tag = path.stem
        df  = pd.read_csv(path, sep="\t")

        # 強制關鍵欄位轉數值
        for col in ["log2FoldChange", "lfcSE", "padj"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        df = df.add_prefix(f"{tag}__")
        df = df.rename(columns={f"{tag}__Geneid": "Geneid"})
        dfs.append(df.set_index("Geneid"))

    if not dfs:
        raise SystemExit(f"[ERROR] No *.tsv files found in {INPUT_DIR}")

    merged = pd.concat(dfs, axis=1, join="inner")

    # ---------- Meta-analysis ----------
    sig_cols = [c for c in merged.columns if c.endswith("__padj")]
    fc_cols  = [c.replace("__padj", "__log2FoldChange") for c in sig_cols]
    se_cols  = [c.replace("__padj", "__lfcSE")          for c in sig_cols]

    merged["meta_p"]      = merged.apply(partial(fisher_p, sig_cols=sig_cols), axis=1)
    merged["meta_log2FC"] = merged.apply(partial(meta_fc, fc_cols=fc_cols, se_cols=se_cols), axis=1)

    padj_mat = merged[sig_cols].fillna(1)           # 不顯著
    fc_mat   = merged[fc_cols].abs().fillna(0)      # 無變化

    tags = [c.split("__", 1)[0] for c in sig_cols]   # 例如 "SRR12345"
    padj_mat.columns = tags
    fc_mat.columns   = tags

    # ---------- Call DEG / Non-DEG ----------
    sig_bool    = (padj_mat < PADJ_TH) & (fc_mat > FC_TH)
    nonsig_bool = (padj_mat >= PADJ_TH) | (fc_mat <= FC_TH)

    merged["sig_count"] = sig_bool.sum(axis=1)
    merged["sig_prop"]  = merged["sig_count"] / len(sig_cols)

    # DEG
    deg = (merged[["sig_count", "sig_prop", "meta_p", "meta_log2FC"]]
           .loc[merged["sig_count"] > 0]
           .dropna(subset=["meta_p"])
           .sort_values(by=["meta_log2FC", "sig_count", "meta_p"],
                        ascending=[False, False, True])
           .reset_index())

    # Non-DEG：所有比較皆非顯著
    mask_non = nonsig_bool.all(axis=1)
    non_deg  = (merged.loc[mask_non,
                           ["sig_count", "sig_prop",
                            "meta_p", "meta_log2FC"]]
                .reset_index())

    # ---------- Save ----------
    TSV_DIR.mkdir(parents=True, exist_ok=True)
    deg.to_csv(OUT_DEG,  sep="\t", index=False)
    non_deg.to_csv(OUT_NON, sep="\t", index=False)

    print(f"[INFO] DEG     → {OUT_DEG}  (n={len(deg):,})")
    print(f"[INFO] Non-DEG → {OUT_NON} (n={len(non_deg):,})")


# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()