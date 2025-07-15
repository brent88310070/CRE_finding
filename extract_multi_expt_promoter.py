#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==========================

● 對 gene list 資料夾樹進行批次 promoter 擷取  
● 每個 <sample>/ 資料夾須至少包含 DEG.txt 與/或 nonDEG.txt  
● 產出：
   <sample>/<sample>_DEG_promoter_1kb.fa
   <sample>/<sample>_nonDEG_promoter_1kb.fa
   <sample>/<sample>_missing_ids.txt  (如有缺基因)

只需修改 CONFIG 區塊即可。
"""

from pathlib import Path
import random, sys, argparse
import pandas as pd
from pyfaidx import Fasta
from tqdm import tqdm

# ---------- 共用工具 ----------
def load_gene_table(gff_path):
    keep_types = {
        "gene", "transposable_element_gene",
        "ncRNA_gene", "lncRNA_gene", "pseudogene"
    }
    rows = []
    with open(gff_path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            seqid, src, ftype, start, end, score, strand, phase, attrs = ln.rstrip().split("\t")
            if ftype not in keep_types:
                continue
            attr = {kv.split("=", 1)[0]: kv.split("=", 1)[1] for kv in attrs.split(";") if "=" in kv}
            gid = attr.get("ID") or attr.get("gene_id") or attr.get("Name")
            if gid:
                rows.append([gid, seqid, int(start), int(end), strand])
    return pd.DataFrame(rows, columns=["gene_id", "chr", "start", "end", "strand"]).set_index("gene_id")


def calc_promoter(row, up_len, chr_len):
    if row.strand == "+":
        p_start = max(1, row.start - up_len)
        p_end = row.start - 1
    else:
        p_start = row.end + 1
        p_end = min(row.end + up_len, chr_len)
    return p_start, p_end


def rc(seq: str) -> str:
    return seq[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))


def write_promoters(id_list, label, sample_name, anno, fa, out_dir):
    if not id_list:
        return 0
    fa_name = out_dir / f"{sample_name}_{label}_promoter_{PROMOTER_UP_BP//1000}kb.fa"
    missing_file = out_dir / f"{sample_name}_{label}_missing_ids.txt"

    with open(fa_name, "w") as fasta_out:
        missing = []
        for gid in tqdm(id_list, desc=f"{sample_name}-{label}", disable=not VERBOSE):
            if gid not in anno.index:
                missing.append(gid)
                continue
            row = anno.loc[gid]
            chr_len = len(fa[row.chr])
            p_start, p_end = calc_promoter(row, PROMOTER_UP_BP, chr_len)
            seq = fa[row.chr][p_start - 1:p_end].seq
            if row.strand == "-":
                seq = rc(seq)
            header = f">{gid}|{row.chr}:{p_start}-{p_end}({row.strand})|{PROMOTER_UP_BP}bp_upstream|{label}"
            fasta_out.write(f"{header}\n{seq}\n")

    if missing:
        missing_file.write_text("\n".join(missing))
        sys.stderr.write(f"[WARN] {sample_name}-{label}: {len(missing)} IDs not in GFF – see {missing_file}\n")
    return len(id_list) - len(missing)


# ---------- 主流程 ----------
def main(argv=None):
    p = argparse.ArgumentParser(description="Extract promoters for many experiments.")
    p.add_argument("-r", "--root_dir",        default="./multi_exp_tomato")
    p.add_argument("--deg_filename",          default="DEG.txt")
    p.add_argument("--nondeg_filename",       default="nonDEG.txt")
    p.add_argument("-g", "--gff",             required=True,
                   help="Genome annotation (GFF/GTF)")
    p.add_argument("-f", "--fasta",           required=True,
                   help="Genome FASTA")
    p.add_argument("-u", "--upstream_bp",     type=int, default=1000,
                   help="Promoter length upstream (bp)")
    p.add_argument("--random_seed",           type=int, default=42)
    p.add_argument("--verbose",               action="store_true")
    args = p.parse_args(argv)

    # 讓 write_promoters 看得到使用者的設定
    global PROMOTER_UP_BP, VERBOSE
    PROMOTER_UP_BP = args.upstream_bp
    VERBOSE        = args.verbose

    if args.random_seed is not None:
        random.seed(args.random_seed)

    root = Path(args.root_dir).expanduser()
    anno = load_gene_table(args.gff)
    fa   = Fasta(args.fasta, sequence_always_upper=True)

    for subdir in sorted(p for p in root.iterdir() if p.is_dir()):
        sample = subdir.name
        deg_path  = subdir / args.deg_filename
        ndeg_path = subdir / args.nondeg_filename
        if not deg_path.exists() and not ndeg_path.exists():
            continue

        deg_ids  = deg_path.read_text().strip().splitlines()  if deg_path.exists()  else []
        ndeg_ids = ndeg_path.read_text().strip().splitlines() if ndeg_path.exists() else []

        n_deg  = write_promoters(deg_ids,  "DEG",    sample, anno, fa, subdir)
        n_ndeg = write_promoters(ndeg_ids, "nonDEG", sample, anno, fa, subdir)

        sys.stderr.write(f"[DONE] {sample}: DEG {n_deg}, non-DEG {n_ndeg}\n")

if __name__ == "__main__":
    main()