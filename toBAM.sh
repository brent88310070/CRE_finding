#!/usr/bin/env bash
# ==== 可調整參數 ==========================================
INPUT_PREFIX="SRP399644_tomato"                 # 專案/實驗前綴
TISSUE="root"                                   # 子資料夾名稱，可留空
BASE_DIR="raw_data"                             # 儲存主資料夾
GENOME_DIR="./ref/S_lycopersicum/star_index"    # STAR 索引路徑
THREADS=12                                      # STAR 執行緒數
PAIRED=false                                     # true=雙端 (_1/_2)；false=單端
GZIP_FASTQ=false                                # 若為 *.fastq.gz 請設 true
# =========================================================

INPUT_DIR="${BASE_DIR}/${INPUT_PREFIX}_${TISSUE}"
OUTPUT_DIR="${BASE_DIR}/${INPUT_PREFIX}_${TISSUE}_bam"
mkdir -p "$OUTPUT_DIR"

# ---------- gzip 讀取參數 ----------
if $GZIP_FASTQ; then
  READ_CMD="--readFilesCommand zcat"
  EXT=".fastq.gz"
else
  READ_CMD=""
  EXT=".fastq"
fi

echo "▶  Alignment mode: $([ "$PAIRED" = true ] && echo 'paired-end' || echo 'single-end')"
echo "   Input dir  : $INPUT_DIR"
echo "   Output dir : $OUTPUT_DIR"
echo "   Genome dir : $GENOME_DIR"
echo "   Threads    : $THREADS"
echo "-------------------------------------------------------"

# ---------- 主流程 ----------
if $PAIRED; then
  #  Paired-end: 迴圈處理 _1.fastq/_2.fastq
  for fq1 in "${INPUT_DIR}"/*_1${EXT}; do
    [[ -e "$fq1" ]] || { echo "找不到 *_1${EXT}；結束。"; break; }
    SAMPLE=$(basename "$fq1" "_1${EXT}")
    fq2="${INPUT_DIR}/${SAMPLE}_2${EXT}"

    if [[ ! -f "$fq2" ]]; then
      echo "⚠  配對檔 $fq2 不存在，跳過 $SAMPLE"
      continue
    fi

    echo "▶  STAR  (paired)  : $SAMPLE"
    STAR --genomeDir "$GENOME_DIR" \
         --readFilesIn "$fq1" "$fq2" \
         $READ_CMD \
         --runThreadN "$THREADS" \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix "${OUTPUT_DIR}/${SAMPLE}_"
    echo "Finished        : $SAMPLE"
  done
else
  #  Single-end: 逐一處理 .fastq
  for fq in "${INPUT_DIR}"/*${EXT}; do
    [[ -e "$fq" ]] || { echo "找不到 *${EXT}；結束。"; break; }
    SAMPLE=$(basename "$fq" "${EXT}")

    echo "▶  STAR  (single) : $SAMPLE"
    STAR --genomeDir "$GENOME_DIR" \
         --readFilesIn "$fq" \
         $READ_CMD \
         --runThreadN "$THREADS" \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix "${OUTPUT_DIR}/${SAMPLE}_"
    echo "Finished        : $SAMPLE"
  done
fi

echo "All samples processed!"
