#!/bin/bash

# ==== 可調整參數 ====
SRR_LIST="srr_list.txt"                  # SRR ID 清單檔案
BASE_DIR="raw_data"                      # 基礎儲存資料夾
THREADS=8                                # fasterq-dump 使用的執行緒數
KEEP_SRA=false                           # 若要保留 .sra 檔，改成 true
# ===================

# 讀取第一行作為資料夾名稱
read PROJECT_NAME < srr_list.txt

# 建立目錄
OUTDIR="raw_data/$PROJECT_NAME"
mkdir -p "$OUTDIR"

# 從第二行開始讀取 SRR ID 並下載
tail -n +2 srr_list.txt | while read SRR; do
    echo "Downloading $SRR ..."
    prefetch "$SRR"
    fasterq-dump "$SRR" --split-files --threads 8 -O "$OUTDIR"
    rm -f "$SRR/$SRR.sra"
    rm -rf "$SRR"
done
