#!/usr/bin/env bash
set -euo pipefail

################################################################
# 🔧 可調整參數區 －－ 只改這裡就能控制整條流程 ###################
################################################################
# ── Step-1：合併各樣本 DEG 結果 ───────────────────────────────
INPUT_DIR="./tomato_deg_results"   # 存放多個 *.tsv
PADJ_TH=0.05                       # FDR 門檻
FC_TH=1                            # |log2FC| 門檻

# ── Step-2：篩選嚴謹 DEG / Non-DEG ───────────────────────────
SPECIES="tomato"                   # 也會作為檔名前綴
SIG_COUNT=3                        # sig_count ≥ ?
DEG_FC_TH=1                        # meta_log2FC > ?
NON_P_TH=0.1                       # Non-DEG 的 meta_p > ?
NON_FC_NEAR0_TH=0.1                # |meta_log2FC| ≤ ?

# ── Step-3：擷取啟動子序列 ────────────────────────────────────
PROMOTER_UP_BP=1000                # 1000 或 2000
NEG_MULTIPLIER=3     # 每個 positive 配幾個 negative
NEG_MIN=1000         # 至少多少條 negative
GFF_PATH="ref/S_lycopersicum/ITAG4.1_gene_models.gff"
FASTA_PATH="ref/S_lycopersicum/S_lycopersicum_chromosomes.4.00.fa"
# GFF_PATH="ref/Araport/Araport11_GFF3_genes_transposons.current.gff"
# FASTA_PATH="ref/Araport/TAIR10_chr_all.fas"
SEED=42                            # Down-sampling Non-DEG 用
################################################################


# --- pick an available python interpreter ---
PYTHON="$(command -v python3 || command -v python)"
if [ -z "$PYTHON" ]; then
  echo "[ERROR] Python not found. Install python3 or activate your conda env first." >&2
  exit 1
fi

echo "=== Step 1. DEG summary ==="
"$PYTHON" deg_summary.py \
  --input_dir "$INPUT_DIR" \
  --padj_th   "$PADJ_TH" \
  --fc_th     "$FC_TH"

echo "=== Step 2. Filter DEG / Non-DEG lists ==="
"$PYTHON" extract_DEG_and_nonDEG.py \
  --species     "$SPECIES" \
  --sig_th      "$SIG_COUNT" \
  --deg_fc_th   "$DEG_FC_TH" \
  --non_p_th    "$NON_P_TH" \
  --non_fc_th   "$NON_FC_NEAR0_TH"

echo "=== Step 3. Extract promoter FASTA ==="
"$PYTHON" extract_promoter.py \
  --sig_count     "$SIG_COUNT" \
  --up_bp         "$PROMOTER_UP_BP" \
  --gff_path      "$GFF_PATH" \
  --fasta_path    "$FASTA_PATH" \
  --prefix        "$SPECIES" \
  --seed          "$SEED" \
  --neg_multiplier "$NEG_MULTIPLIER" \
  --neg_min        "$NEG_MIN"

echo "✅  Pipeline finished – promoter FASTA files are in ./prom_seq_files/"
