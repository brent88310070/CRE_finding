#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
filter_deg_lists.py
===================
Apply stringent thresholds to merged DEG / Non-DEG tables.

輸入檔：
  ./deg_summary/{species}_DEG.tsv
  ./deg_summary/{species}_nonDEG.tsv
輸出檔：
  ./deg_summary/{species}_DEG_filtered_sig_count_{sig_th}.tsv
  ./deg_summary/{species}_nonDEG_filtered_sig_count_{sig_th}.tsv
以及對應 _geneid.txt 清單
"""

from pathlib import Path
import argparse
import pandas as pd


# ----------------------------------------------------------------------
# Argument parser
# ----------------------------------------------------------------------
def get_args():
    p = argparse.ArgumentParser(
        description="Filter rigorous DEG / Non-DEG gene lists"
    )
    p.add_argument("--species",     default="tomato",
                   help="File prefix of species (default: tomato)")
    p.add_argument("--summary_dir", default="./deg_summary",
                   help="Directory containing merged summary TSVs")
    p.add_argument("--sig_th",      type=int,   default=3,
                   help="Min. sig_count for DEG")
    p.add_argument("--deg_fc_th",   type=float, default=1.0,
                   help="meta_log2FC > this for DEG")
    p.add_argument("--non_p_th",    type=float, default=0.1,
                   help="meta_p > this for Non-DEG")
    p.add_argument("--non_fc_th",   type=float, default=0.1,
                   help="|meta_log2FC| ≤ this for Non-DEG")
    return p.parse_args()


# ----------------------------------------------------------------------
def main():
    args = get_args()

    sp          = args.species
    SUM_DIR     = Path(args.summary_dir)
    IN_DEG      = SUM_DIR / f"{sp}_DEG.tsv"
    IN_NON      = SUM_DIR / f"{sp}_nonDEG.tsv"

    # 門檻
    SIG_TH      = args.sig_th
    DEG_FC_TH   = args.deg_fc_th
    NON_P_TH    = args.non_p_th
    NON_FC_TH   = args.non_fc_th

    OUT_DEG = SUM_DIR / f"{sp}_DEG_filtered_sig_count_{SIG_TH}.tsv"
    OUT_NON = SUM_DIR / f"{sp}_nonDEG_filtered_sig_count_{SIG_TH}.tsv"

    # ---------- 1. 讀檔 ----------
    deg  = pd.read_csv(IN_DEG, sep="\t")
    non  = pd.read_csv(IN_NON, sep="\t")

    # ---------- 2. 過濾 ----------
    # 2-1 DEG：sig_count ≥ SIG_TH 且 meta_log2FC > DEG_FC_TH
    deg_filt = (deg
        .loc[(deg["sig_count"] >= SIG_TH) & (deg["meta_log2FC"] > DEG_FC_TH)]
        .reset_index(drop=True))

    # 2-2 Non-DEG：meta_p > NON_P_TH 且 |meta_log2FC| ≤ NON_FC_TH
    non_filt = (non
        .dropna(subset=["meta_p", "meta_log2FC"])
        .loc[(non["meta_p"] > NON_P_TH) &
             (non["meta_log2FC"].abs() <= NON_FC_TH)]
        .reset_index(drop=True))

    # ---------- 3. 輸出 ----------
    deg_filt.to_csv(OUT_DEG,  sep="\t", index=False)
    non_filt.to_csv(OUT_NON, sep="\t", index=False)

    # ---------- 4. 摘要 ----------
    print(f"[INFO] DEG  檔案：{OUT_DEG}   (n={len(deg_filt):,})")
    print(deg_filt.head(), "\n")

    print(f"[INFO] Non-DEG 檔案：{OUT_NON} (n={len(non_filt):,})")
    print(non_filt.head())

    # ---------- 5. 輸出純 Geneid 清單 ----------
    deg_filt["Geneid"].to_csv(
        OUT_DEG.with_name(OUT_DEG.stem + "_geneid.txt"),
        index=False, header=False
    )
    non_filt["Geneid"].to_csv(
        OUT_NON.with_name(OUT_NON.stem + "_geneid.txt"),
        index=False, header=False
    )


# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()
