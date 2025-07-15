#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch split DEG / non-DEG and save Geneid lists.
"""

from pathlib import Path
import pandas as pd
import random
import argparse
import sys


# ────── 主工具函式 ──────
def split_deg(
    file_path: Path,
    out_root: Path,
    basemean_deg_min: float,
    log2fc_deg_min: float,
    padj_deg_max: float,
    log2fc_non_max: float,
    padj_non_min: float,
    neg_multiplier: int,
    neg_min: int,
    random_seed: int | None,
    verbose: bool,
):
    """將單一 DESeq2 TSV 拆分成 DEG 與 non-DEG Geneid 名單"""
    df = pd.read_csv(file_path, sep="\t")

    # ------------ 篩選正樣本 (DEG) ------------
    deg_mask = (
        (df["baseMean"] >= basemean_deg_min)
        & (df["log2FoldChange"] > log2fc_deg_min)
        & (df["padj"] < padj_deg_max)
    )
    deg_ids = df.loc[deg_mask, "Geneid"].dropna().astype(str)

    # ------------ 篩選負樣本 (non-DEG) ------------
    non_mask = (df["log2FoldChange"] < log2fc_non_max) & (
        df["padj"] > padj_non_min
    )
    non_ids_all = df.loc[non_mask, "Geneid"].dropna().astype(str)

    # ------------ 抽樣負樣本：3× 或至少 1000 ------------
    if random_seed is not None:
        random.seed(random_seed)

    target_neg = max(len(deg_ids) * neg_multiplier, neg_min)
    sample_size = min(len(non_ids_all), target_neg)
    non_ids = pd.Series(random.sample(list(non_ids_all), sample_size))

    # ------------ 輸出 ------------
    out_dir = out_root / file_path.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    deg_ids.to_csv(out_dir / "DEG.txt", index=False, header=False)
    non_ids.to_csv(out_dir / "nonDEG.txt", index=False, header=False)

    # ------------ 訊息 ------------
    if verbose:
        print(
            f"[✓] {file_path.name:>20}  "
            f"DEG={len(deg_ids):5d} | "
            f"non-DEG={len(non_ids):5d} "
            f"(requested {target_neg}, original {len(non_ids_all)})"
        )


# ────── 主程式 ──────
def main(argv=None):
    p = argparse.ArgumentParser(
        description="Split each DESeq2 TSV into DEG / non-DEG Geneid lists."
    )
    p.add_argument("-i", "--input_dir", default="./tomato_deg_results")
    p.add_argument("-o", "--output_dir", default="./multi_exp_tomato")
    p.add_argument("--file_pattern", default="*.tsv")

    # DEG 門檻
    p.add_argument("--basemean_deg_min", type=float, default=10)
    p.add_argument("--log2fc_deg_min", type=float, default=1)
    p.add_argument("--padj_deg_max", type=float, default=0.05)

    # non-DEG 門檻
    p.add_argument("--log2fc_non_max", type=float, default=0.1)
    p.add_argument("--padj_non_min", type=float, default=0.1)

    # 負樣本數量規則
    p.add_argument(
        "--neg_multiplier",type=int,default=3,
        help="Negatives per positive (default: 3)",)
    p.add_argument(
        "--neg_min",type=int,default=1000,
        help="Minimum number of negatives to keep (default: 1000)",)

    p.add_argument("--random_seed", type=int, default=42)
    p.add_argument("--verbose", action="store_true")

    args = p.parse_args(argv)

    in_root = Path(args.input_dir).expanduser().resolve()
    out_root = Path(args.output_dir).expanduser().resolve()
    out_root.mkdir(exist_ok=True)

    tsv_files = sorted(in_root.glob(args.file_pattern))
    if not tsv_files:
        sys.exit(f"No files matched '{args.file_pattern}' in {in_root}")

    for fp in tsv_files:
        split_deg(
            fp,
            out_root,
            args.basemean_deg_min,
            args.log2fc_deg_min,
            args.padj_deg_max,
            args.log2fc_non_max,
            args.padj_non_min,
            args.neg_multiplier,
            args.neg_min,
            args.random_seed,
            args.verbose,
        )


if __name__ == "__main__":
    main()
