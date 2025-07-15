#!/usr/bin/env bash
# ========================================
POS=./prom_seq_files/tomato_DEG_promoter_1kb_sig_count_3.fa          # 正集合
NEG=./prom_seq_files/tomato_nonDEG_promoter_1kb_sig_count_3.fa       # 負集合
OUT=motif_out
OUT_DIR=tomato_1kb_sig_count_3
MINW=5            # 最短 motif 長度
MAXW=15           # 最長 motif 長度
N_MOTIFS=20
# ========================================

mkdir -p "$OUT"
mkdir -p "$OUT/$OUT_DIR"

echo "--- Run STREME ---"
streme \
  --p  "$POS" \
  --n  "$NEG" \
  --dna \
  --minw $MINW --maxw $MAXW \
  --nmotifs $N_MOTIFS \
  --oc "$OUT/$OUT_DIR" \
  --verbosity 1

echo "✅ Motif analysis finished. See $OUT/index.html"