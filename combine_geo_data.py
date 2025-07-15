import pandas as pd
import re

# ========== 可調整參數區 ==========
file_name = './raw_data/SRP399644_tomato_root_count/SRP399644_tomato_counts.txt'    # featureCounts 輸出
output_name = './exp_files/tomato/SRP399644_tomato_root_exp.tsv'                    # 最終輸出檔案
meta_path = 'run_info.txt'                                                          # 需包含 SRR 對 GEO ID 的對照表
# ==================================

# 讀取 featureCounts 結果，略過 # 開頭的註解行
df = pd.read_csv(file_name, sep="\t", comment="#")
df = df.set_index("Geneid")

# 只保留 counts 欄位（第六欄以後）
cts = df.iloc[:, 5:]

# 從欄位名稱中擷取 SRR 編號
cts.columns = [
    re.search(r'SRR\d+', col).group(0) if re.search(r'SRR\d+', col) else col
    for col in cts.columns
]

# 讀取 SRR 與 GEO 對照表，並建立 ID 對應關係
meta = pd.read_csv(meta_path, sep="\t")
id_map = dict(zip(meta["Run"], meta["GEO_Accession (exp)"]))

# 替換欄名（SRR → GEO ID），並合併同一 GEO ID 的樣本
cts_renamed = cts.rename(columns=id_map)
cts_combined = cts_renamed.T.groupby(level=0).sum().T


# 輸出為 TSV
cts_combined.to_csv(output_name, sep="\t")
