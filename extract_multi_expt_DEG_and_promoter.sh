#!/usr/bin/env bash
# =====================================================================

#   1) 依 DESeq2 結果產生每個實驗的 DEG / non-DEG GeneID 清單
#   2) 批次擷取各實驗 –k bp 範圍的 promoter FASTA
# ---------------------------------------------------------------------
# 使用前請確認：
#   • extract_multi_expt_DEG_and_nonDEG.py
#   • extract_multi_expt_promoter.py
#   兩支腳本已放在同一資料夾
# =====================================================================

set -euo pipefail         # 任一步驟出錯即終止腳本

# ────────────────────── CONFIG 可調參數 ──────────────────────
# 1) DEG / non-DEG GeneID 清單 (extract_multi_expt_DEG_and_nonDEG.py)
INPUT_DIR="./ara_deg_results"          # 存放 DESeq2 .tsv 的目錄
OUTPUT_DIR="./multi_exp_arabidopsis"           # 產出清單的根目錄
FILE_PATTERN="*.tsv"                      # 要讀取的檔案格式

BASEMEAN_DEG_MIN=10                       # DEG 條件
LOG2FC_DEG_MIN=1
PADJ_DEG_MAX=0.05

LOG2FC_NON_MAX=0.1                        # non-DEG 條件
PADJ_NON_MIN=0.1
NEG_MULTIPLIER=3     # 每個 positive 配幾個 negative
NEG_MIN=1000         # 至少多少條 negative

# 2) Promoter 擷取 (extract_multi_expt_promoter.py)
GFF_PATH="ref/Araport/Araport11_GFF3_genes_transposons.current.gff"
FASTA_PATH="ref/Araport/TAIR10_chr_all.fas"
PROMOTER_UP_BP=1000                       # 1000 = −1 kb，2000 = −2 kb …
DEG_FILENAME="DEG.txt"                    # 若改名請同步修改
NONDEG_FILENAME="nonDEG.txt"

# 3) 其他
RANDOM_SEED=42
VERBOSE=true                              # 是否顯示進度
# ───────────────────────────────────────────────────────────

# ===== (1) 產生 DEG / non-DEG GeneID 清單 =====
echo "=== Step 1. Get DEG / non-DEG GeneID list ==="
python extract_multi_expt_DEG_and_nonDEG.py \
    --input_dir          "${INPUT_DIR}" \
    --output_dir         "${OUTPUT_DIR}" \
    --file_pattern       "${FILE_PATTERN}" \
    --basemean_deg_min   "${BASEMEAN_DEG_MIN}" \
    --log2fc_deg_min     "${LOG2FC_DEG_MIN}" \
    --padj_deg_max       "${PADJ_DEG_MAX}" \
    --log2fc_non_max     "${LOG2FC_NON_MAX}" \
    --padj_non_min       "${PADJ_NON_MIN}" \
    --neg_multiplier     "${NEG_MULTIPLIER}" \
    --neg_min            "${NEG_MIN}" \
    --random_seed        "${RANDOM_SEED}" \
    $( [[ "${VERBOSE}" == true ]] && echo "--verbose" )

# ===== (2) 批次擷取 promoter FASTA =====
echo "=== Step 2. Extract promoter FASTA ==="
python extract_multi_expt_promoter.py \
    --root_dir        "${OUTPUT_DIR}" \
    --deg_filename    "${DEG_FILENAME}" \
    --nondeg_filename "${NONDEG_FILENAME}" \
    --gff             "${GFF_PATH}" \
    --fasta           "${FASTA_PATH}" \
    --upstream_bp     "${PROMOTER_UP_BP}" \
    --random_seed     "${RANDOM_SEED}" \
    $( [[ "${VERBOSE}" == true ]] && echo "--verbose" )

echo -e "\n Finished successfully."
