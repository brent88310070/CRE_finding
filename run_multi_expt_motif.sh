#!/usr/bin/env bash
#========================================================================
# run_streme_batch.sh – batch motif search for every sample folder
#
# 需先安裝 MEME Suite 並將 `streme` 加入 $PATH
# 執行：  bash run_streme_batch.sh
#========================================================================

######################### CONFIG (edit here) ############################
ROOT_DIR="./multi_exp_arabidopsis"            # ← 放各 sample 子資料夾的根目錄
PROM_SIZE="1kb"                  # "1kb" / "2kb"... (檔名必須含此字串)
MINW=5                           # 最短 motif 長度
MAXW=15                          # 最長 motif 長度
N_MOTIFS=20                      # 要找幾個 motif
VERBOSITY=1                      # streme --verbosity
#########################################################################

echo -e "\n=== STREME batch run ==="
echo   "Root dir      : $ROOT_DIR"
echo   "Promoter size : $PROM_SIZE"
echo   "Motif length  : $MINW–$MAXW"
echo   "Motif number  : $N_MOTIFS"
echo

shopt -s nullglob

for sample_dir in "$ROOT_DIR"/SRP*/ ; do
    sample=$(basename "$sample_dir")

    # 自動尋找符合命名規則的 FASTA；用變數而非陣列更容易 debug
    pos_file=$(echo "$sample_dir"/*_DEG_promoter_"$PROM_SIZE".fa | head -n1)
    neg_file=$(echo "$sample_dir"/*_nonDEG_promoter_"$PROM_SIZE".fa | head -n1)

    if [[ -f "$pos_file" && -f "$neg_file" ]]; then
        oc_dir="$sample_dir/streme_${PROM_SIZE}"
        mkdir -p "$oc_dir"

        echo "--- [$sample] Run STREME ---"
        streme \
          --p "$pos_file" \
          --n "$neg_file" \
          --dna \
          --minw "$MINW" --maxw "$MAXW" \
          --nmotifs "$N_MOTIFS" \
          --oc "$oc_dir" \
          --verbosity "$VERBOSITY"

        echo "✅ [$sample] Motif search done → $oc_dir"
        echo
    else
        echo "⚠️  [$sample] FASTA not found – skipped."
    fi
done

echo "=== All done ==="
