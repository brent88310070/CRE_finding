#!/usr/bin/env python3
# ================================================================
# merge_and_deseq2.py
# 1. 將 ./exp_files/tomato/*.tsv 合併成 Tomato_all_exp.tsv
# 2. 依 read_treat_control_list.txt 指定的配對，批次執行 DESeq2
# ------------------------------------------------
# 需要 pip install pydeseq2 pandas
# ================================================================
import pandas as pd
import glob, os, pathlib
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# ==========================
# 可調整參數區
# ==========================
# --- 合併 expression 檔 ---
input_folder   = "./exp_files/tomato"        # 單實驗 expression 檔資料夾
file_pattern   = "*.tsv"                     # 檔案格式
index_col      = "Geneid"                    # 基因欄位名稱
merged_counts  = "./exp_files/Tomato_all_exp.tsv"

# --- DESeq2 相關 ---
sample_info_path = "./read_treat_control_list.txt"  # Treatment/Control 配對設定
deg_output_dir   = "./tomato_deg_results"           # DESeq2 輸出資料夾
# ===============================================================


# ---------- Step 1：合併 expression 表 ----------
print("► 合併 expression 檔 ...")
file_list = glob.glob(os.path.join(input_folder, file_pattern))
if not file_list:
    raise FileNotFoundError(f"No files matched {file_pattern} in {input_folder}")

merged_df = None
for i, f in enumerate(sorted(file_list)):
    df = pd.read_csv(f, sep="\t").set_index(index_col)
    merged_df = df if merged_df is None else merged_df.join(df, how="outer")
    print(f"  ({i+1:>2}/{len(file_list)}) {os.path.basename(f)} merged")

merged_df.to_csv(merged_counts, sep="\t")
print(f"✔ 合併完成：{merged_counts}  (genes={merged_df.shape[0]}, samples={merged_df.shape[1]})\n")


# ---------- Step 2：解析 Treatment/Control 配對 ----------
def parse_txt_to_list_of_dict(file_path):
    """讀取兩行一組的配對設定檔，回傳 list[dict]."""
    df = pd.read_csv(file_path, sep="\t", header=None).dropna(how="all").reset_index(drop=True)
    result = []
    for i in range(0, len(df), 2):
        dataset_id = df.iloc[i, 0]
        control    = df.iloc[i,   2:].dropna().tolist()
        treatment  = df.iloc[i+1, 2:].dropna().tolist()
        result.append({"dataset_id": dataset_id, "control": control, "treatment": treatment})
    return result


dataset_list = parse_txt_to_list_of_dict(sample_info_path)
print(f"► 共有 {len(dataset_list)} 個資料集準備進行 DESeq2\n")


# ---------- Step 3：批次跑 PyDESeq2 ----------
counts_df = pd.read_csv(merged_counts, sep="\t", index_col=0)
deg_dir   = pathlib.Path(deg_output_dir)
deg_dir.mkdir(parents=True, exist_ok=True)

for n, entry in enumerate(dataset_list, 1):
    ds_id     = entry["dataset_id"]
    control   = entry["control"]
    treatment = entry["treatment"]

    print(f"[{n}/{len(dataset_list)}] {ds_id}: {len(treatment)} treat vs {len(control)} ctrl")

    meta = (
        pd.DataFrame({
            "sample":    treatment + control,
            "condition": ["treatment"] * len(treatment) + ["control"] * len(control)
        })
        .set_index("sample")
    )

    # 建立 DESeq2 dataset
    dds = DeseqDataSet(
        counts=counts_df.loc[:, meta.index].T,   # DESeq2 需要 sample × gene
        metadata=meta,
        design_factors="condition",
        ref_level=["condition", "control"]
    )
    dds.deseq2()

    stats = DeseqStats(dds, contrast=["condition", "treatment", "control"])
    stats.summary()

    deg = (
        stats.results_df
        .sort_values(["padj", "log2FoldChange"], ascending=[True, False])
        .reset_index()
        .rename(columns={"index": "Gene"})
    )

    out_path = deg_dir / f"{ds_id}.tsv"
    deg.to_csv(out_path, sep="\t", index=False)
    print(f"  → DEGs saved: {out_path}")

print(f"\n 全部完成！結果輸出於：{deg_output_dir}")
