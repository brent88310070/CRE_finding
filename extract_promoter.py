#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_promoter.py
===================
Grab −1 kb / −2 kb promoter sequences for DEG & Non-DEG gene lists.

輸入：
  ./deg_summary/{prefix}_DEG_filtered_sig_count_{sig_count}_geneid.txt
  ./deg_summary/{prefix}_nonDEG_filtered_sig_count_{sig_count}_geneid.txt
"""

import sys, random, argparse
from pathlib import Path
import pandas as pd
from pyfaidx import Fasta
from tqdm import tqdm


# ----------------------------------------------------------------------
# Argument parser
# ----------------------------------------------------------------------
def get_args():
    p = argparse.ArgumentParser(
        description="Extract upstream promoter sequences for DEG / Non-DEG"
    )
    p.add_argument("--sig_count", type=int,   default=3,
                   help="sig_count threshold used in previous step")
    p.add_argument("--neg_multiplier", type=int, default=3,
               help="How many negatives per positive (default: 3)")
    p.add_argument("--neg_min", type=int, default=1000,
                help="Minimum number of negatives to keep (default: 1000)")
    p.add_argument("--summary_dir", default="deg_summary",
                   help="Folder containing *_geneid.txt produced in Step-2")
    p.add_argument("--gff_path", default="ref/S_lycopersicum/ITAG4.1_gene_models.gff",
                   help="GFF3 / GTF annotation file")
    p.add_argument("--fasta_path", default="ref/S_lycopersicum/S_lycopersicum_chromosomes.4.00.fa",
                   help="Reference genome FASTA")
    p.add_argument("--up_bp", type=int, default=1000,
                   help="Upstream length (bp) to extract (e.g. 1000 or 2000)")
    p.add_argument("--prefix",  default="tomato",
                   help="Species / file prefix")
    p.add_argument("--seed",    type=int, default=42,
                   help="Random seed for down-sampling Non-DEG (None = off)")
    p.add_argument("--out_dir", default="prom_seq_files",
                   help="Output folder for FASTA & log files")
    return p.parse_args()


# ----------------------------------------------------------------------
def read_id_list(path):
    if not Path(path).exists():
        return []
    return [ln.strip() for ln in open(path) if ln.strip()]


def load_gene_table(gff_path):
    keep_types = {
        "gene", "transposable_element_gene",
        "ncRNA_gene", "lncRNA_gene", "pseudogene",
    }
    rows = []
    with open(gff_path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            seqid, src, ftype, start, end, score, strand, phase, attrs = ln.rstrip().split("\t")
            if ftype not in keep_types:
                continue
            attr = {kv.split("=", 1)[0]: kv.split("=", 1)[1]
                    for kv in attrs.split(";") if "=" in kv}
            gid = attr.get("ID") or attr.get("gene_id") or attr.get("Name")
            if gid is None:
                continue
            rows.append([gid, seqid, int(start), int(end), strand])
    df = pd.DataFrame(rows, columns=["gene_id", "chr", "start", "end", "strand"])
    return df.set_index("gene_id")


def calc_promoter(row, up_len, chr_len):
    if row.strand == "+":
        p_start = max(1, row.start - up_len)
        p_end = row.start - 1
    else:
        p_start = row.end + 1
        p_end = min(row.end + up_len, chr_len)
    return p_start, p_end


def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))


# ----------------------------------------------------------------------
def main():
    args = get_args()

    SIG_COUNT      = args.sig_count
    NEG_MULTIPLIER = args.neg_multiplier
    NEG_MIN        = args.neg_min
    PROMOTER_UP_BP = args.up_bp
    PREFIX         = args.prefix
    SEED           = args.seed
    OUT_DIR        = args.out_dir

    Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

    # Input gene-id list paths
    DEG_FILE    = f"{args.summary_dir}/{PREFIX}_DEG_filtered_sig_count_{SIG_COUNT}_geneid.txt"
    NONDEG_FILE = f"{args.summary_dir}/{PREFIX}_nonDEG_filtered_sig_count_{SIG_COUNT}_geneid.txt"

    # Reproducible down-sampling
    if SEED is not None:
        random.seed(SEED)

    # Read gene lists
    deg_ids      = read_id_list(DEG_FILE)
    nondeg_all   = read_id_list(NONDEG_FILE)

    # Down-sample Non-DEG to match DEG count
    target_neg = max(len(deg_ids) * NEG_MULTIPLIER, NEG_MIN)

    # ❷ 依照目標數量抽樣 Non-DEG
    if deg_ids and nondeg_all and len(nondeg_all) >= target_neg:
        nondeg_ids = random.sample(nondeg_all, k=target_neg)
        sys.stderr.write(
            f"[INFO] Down-sampled {target_neg} Non-DEG IDs "
            f"(3× positives, ≥{NEG_MIN}) out of {len(nondeg_all)}\n"
        )
    else:
        nondeg_ids = nondeg_all       # 不足目標數，全部保留
        sys.stderr.write(
            f"[INFO] Only {len(nondeg_all)} Non-DEG IDs available; kept them all\n"
        )

    # Load annotation & genome
    anno = load_gene_table(args.gff_path)
    fa   = Fasta(args.fasta_path, sequence_always_upper=True)

    # --------------------------------------------------------------
    def write_promoters(id_list, label):
        if not id_list:
            return 0

        out_path = (Path(OUT_DIR) /
                    f"{PREFIX}_{label}_promoter_{PROMOTER_UP_BP//1000}kb_sig_count_{SIG_COUNT}.fa")
        missing = []

        with out_path.open("w") as out_fa:
            for gid in tqdm(id_list, desc=label):
                if gid not in anno.index:
                    missing.append(gid)
                    continue
                row = anno.loc[gid]
                chr_len = len(fa[row.chr])
                p_start, p_end = calc_promoter(row, PROMOTER_UP_BP, chr_len)
                seq = fa[row.chr][p_start - 1: p_end].seq
                if row.strand == "-":
                    seq = reverse_complement(seq)
                header = (f">{gid}|{row.chr}:{p_start}-{p_end}({row.strand})|"
                          f"{PROMOTER_UP_BP}bp_upstream|{label}")
                out_fa.write(f"{header}\n{seq}\n")

        # Write missing ID log
        if missing:
            miss_path = (Path(OUT_DIR) /
                         f"{PREFIX}_{label}_sig_count_{SIG_COUNT}_missing_ids.txt")
            Path(miss_path).write_text("\n".join(missing))
            sys.stderr.write(
                f"[WARN] {label}: {len(missing)} IDs not in GFF (see {miss_path})\n"
            )
        return len(id_list) - len(missing)

    # --------------------------------------------------------------
    n_deg  = write_promoters(deg_ids, "DEG") if deg_ids else 0
    n_ndeg = write_promoters(nondeg_ids, "nonDEG") if nondeg_ids else 0

    sys.stderr.write(
        f"[DONE] Extracted {n_deg} DEG promoters and {n_ndeg} Non-DEG promoters.\n"
    )


# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()
