#!/usr/bin/env bash

# ==== 可調整參數 ==========================================
INPUT_PREFIX="SRP399644_tomato"                 # 專案名稱
TISSUE="root"                                   # 子資料夾名稱
BASE_DIR="raw_data"                             # 儲存主資料夾
THREADS=12                                      # 使用的執行緒數
GTF="ref/S_lycopersicum/ITAG4.1_gene_models.gtf"     # 註解檔案（GTF）
IS_PAIRED=false                                  # 是否為 paired-end 資料
# =========================================================

INPUT_DIR="${BASE_DIR}/${INPUT_PREFIX}_${TISSUE}_bam"
OUTPUT_DIR="${BASE_DIR}/${INPUT_PREFIX}_${TISSUE}_count"
OUTPUT="${OUTPUT_DIR}/${INPUT_PREFIX}_counts.txt"

mkdir -p "$OUTPUT_DIR"

# ==== 收集 BAM 檔 ====
mapfile -t BAM_FILES < <(find "$INPUT_DIR" -type f -name "*Aligned.sortedByCoord.out.bam" | sort)

echo "找到 ${#BAM_FILES[@]} 個 BAM 檔："
for f in "${BAM_FILES[@]}"; do
    echo "  - $f"
done

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "❌ 沒有找到任何 BAM 檔，請確認 INPUT_DIR 是否正確：$INPUT_DIR" >&2
    exit 1
fi

# ==== 組合 featureCounts 指令 ====
FEATURECOUNTS_CMD=(
    featureCounts
    -T "$THREADS"
    -t exon
    -g gene_id
    -a "$GTF"
    -o "$OUTPUT"
)

# 如果是 paired-end，加上 -p
if [ "$IS_PAIRED" = true ]; then
    FEATURECOUNTS_CMD+=(-p)
fi

# 固定使用 primary 對齊結果
FEATURECOUNTS_CMD+=(--primary)

# ==== 執行 featureCounts ====
echo "執行 featureCounts..."
echo "Started at: $(date)"
time "${FEATURECOUNTS_CMD[@]}" "${BAM_FILES[@]}"
echo "Finished at: $(date)"
echo "完成，輸出：$OUTPUT"