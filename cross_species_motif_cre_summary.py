#!/usr/bin/env python3
# ----------------------------------------------------------
# cross_species_motif_cre_summary.py
# 交叉比對兩物種的 filtered.meme or streme.txt，找出共同 motif
# ----------------------------------------------------------
import os, sys, shutil, subprocess
from collections import defaultdict, deque

# ========= 可調整參數 =========
SUMMARY_TYPE = "deg_summary"  # 或 "cre_summary"

TOMTOM_THRESH   = 0.05                       # Tomtom q-value 上限
TEMP_DIR        = "tomtom_cross_temp"        # Tomtom 暫存資料夾
OUT_DIR         = "cross_species_motif"

if SUMMARY_TYPE == "deg_summary":
    SPECIES1_FILE   = "./motif_out/arabidopsis_1kb_sig_count_6/streme.txt"
    SPECIES2_FILE   = "./motif_out/tomato_1kb_sig_count_1/streme.txt"
    CRE_DIR = os.path.join(OUT_DIR, "deg_summary")
elif SUMMARY_TYPE == "cre_summary":
    SPECIES1_FILE   = "./multi_exp_arabidopsis/filtered.meme"
    SPECIES2_FILE   = "./multi_exp_tomato/filtered.meme"
    CRE_DIR = os.path.join(OUT_DIR, "cre_summary")
OUT_TSV         = os.path.join(CRE_DIR, "repeat_motif_cross_species.tsv")
# ============================


# ---------- 共用工具 ----------

def parse_meme_file(meme_file: str):
    """回傳 (header_lines, {motif_id: full_text}, {motif_id: e-value})"""
    header_lines, motifs, evalues = [], {}, {}
    try:
        with open(meme_file) as fh:
            txt = fh.read()
    except Exception as e:
        sys.exit(f"❌ 讀取 {meme_file} 失敗：{e}")

    first = txt.find("MOTIF")
    header_lines = txt[:first].strip().splitlines() if first != -1 else txt.strip().splitlines()
    for blk in txt.split("MOTIF ")[1:]:
        lines = blk.strip().splitlines()
        motif_full = "MOTIF " + blk
        motif_id   = lines[0].split(' STREME-')[0]  # 去掉 STREME-1 之類
        e_val = float("inf")
        for ln in lines:
            if "E=" in ln:
                try:
                    e_val = float(ln.split("E=")[1].strip())
                    break
                except ValueError:
                    pass
        motifs[motif_id]  = motif_full
        evalues[motif_id] = e_val
    return header_lines, motifs, evalues


def run_tomtom(query: str, target: str, out_dir: str, thresh: float):
    """執行 Tomtom，比對 query vs target，回傳 tsv 路徑"""
    cmd = [
        "tomtom",
        "--oc", out_dir,
        "--thresh", str(thresh),
        "--dist", "pearson",
        query, target
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return os.path.join(out_dir, "tomtom.tsv")


def build_cross_clusters(tomtom_tsv: str, evalues: dict):
    """
    把 tomtom 連線視為無向圖，取連通元件。
    代表 motif = 其 component 中 E-value 最小者 (同分取字典序最小)。
    回傳 rep_to_dups {rep: [other_motifs]}
    """
    graph = defaultdict(set)
    with open(tomtom_tsv) as fh:
        for ln in fh:
            if ln.startswith('#') or not ln.strip():
                continue
            if ln.lower().startswith('query_id'):
                continue
            q, t = ln.split('\t', 2)[:2]
            if q == t:
                continue
            graph[q].add(t)
            graph[t].add(q)

    visited, rep_to_dups = set(), {}
    for node in graph:
        if node in visited:
            continue
        comp = []
        dq = deque([node])
        while dq:
            cur = dq.popleft()
            if cur in visited:
                continue
            visited.add(cur)
            comp.append(cur)
            dq.extend(graph[cur] - visited)

        if len(comp) == 1:           # 只有一個 → 無 cross-match
            continue
        comp.sort(key=lambda m: (evalues.get(m, float('inf')), m))
        rep, dups = comp[0], comp[1:]
        rep_to_dups[rep] = dups
    return rep_to_dups


def write_summary(rep_to_dups: dict, evalues: dict, out_tsv: str):
    """寫出 cross-species motif 重複摘要"""
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    reps = sorted(rep_to_dups, key=lambda m: (evalues.get(m, float('inf')), m))
    with open(out_tsv, "w") as fh:
        fh.write("representative\te_value\tdup_count\tduplicate_motifs\n")
        for rep in reps:
            e_val = evalues.get(rep, float('nan'))
            fh.write(
                f"{rep}\t{e_val:.3g}\t{len(rep_to_dups[rep])}\t{';'.join(rep_to_dups[rep])}\n"
            )
    print(f"➡ 已輸出交叉物種重複 motifs → {out_tsv}")


# ----------------  主程式  ----------------
def main():
    # 0) 準備目錄
    shutil.rmtree(TEMP_DIR, ignore_errors=True)
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1) 解析兩個 MEME
    _, motifs1, evals1 = parse_meme_file(SPECIES1_FILE)
    _, motifs2, evals2 = parse_meme_file(SPECIES2_FILE)
    all_evals = {**evals1, **evals2}   # 合併字典

    # 2) Tomtom cross-comparison
    tsv = run_tomtom(SPECIES1_FILE, SPECIES2_FILE, TEMP_DIR, TOMTOM_THRESH)

    # 3) 組 cross-species clusters
    rep_to_dups = build_cross_clusters(tsv, all_evals)

    # 4) 輸出摘要
    write_summary(rep_to_dups, all_evals, OUT_TSV)

    print(f"✅ 完成！檢測到 {len(rep_to_dups)} 組跨物種共通 motifs")


if __name__ == "__main__":
    main()
