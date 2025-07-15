#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================
(1) 掃描 STREME 結果 → 轉檔、E-value 過濾並重新命名
(2) 合併成 all.meme
(3) Tomtom 自比對去冗餘
(4) 繪製 filtered.meme 前 N_LOGO_MOTIFS 個 logo
(5) 依 redundant_motif_reps.tsv 再畫代表 motif 前 N_REP_LOGOS 個 logo
"""

# ───────────────────────── CONFIG ─────────────────────────
# ── 資料來源與輸出 ──
DATA_DIR        = "./multi_exp_tomato"          # 多個實驗資料夾

# ── 去冗餘後輸出 ──
FILTERED_MEME   = f"{DATA_DIR}/filtered.meme"
KEPT_ID_TSV     = f"{DATA_DIR}/kept_motif_ids.tsv"
REDUNDANT_TSV   = f"{DATA_DIR}/redundant_motif_reps.tsv"

# ── 參數設定 ──
EVALUE_FILTER   = 10.0        # STREME motifs：保留 E ≤ 10
TOMTOM_THRESH   = 0.05        # Tomtom q-value 閾值

# ── Logo 繪圖 ──
N_LOGO_MOTIFS   = 5                       # filtered.meme 前 N 個
LOGO_OUT_DIR    = f"{DATA_DIR}/motif_logos"

N_REP_LOGOS     = 5                       # 代表 motif 前 N
REP_LOGO_DIR    = f"{DATA_DIR}/reps_motif_logos"

FIGSIZE         = (4, 1.5)                # 單圖尺寸 (inch)
DPI             = 200                     # 解析度
COLOR_SCHEME    = "classic"               # Logomaker 色盤

# ── 其他 (通常無須調整) ──
OUT_DIR         = f"{DATA_DIR}/temp"            # 暫存目錄
PREPARED_DIR    = f"{OUT_DIR}/prepared"         # 個別實驗過濾後 .meme
ALL_MEME        = f"{OUT_DIR}/all.meme"         # 合併檔
TEMP_DIR        = f"{OUT_DIR}/tomtom_temp"
# ──────────────────────────────────────────────────────────


import os, re, shutil, subprocess, sys
from collections import defaultdict, deque
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker as lm


# ╭────────────────────── 共用工具 ───────────────────────╮
def run(cmd, **kw):
    try:
        return subprocess.run(cmd, check=True, text=True,
                              capture_output=True, **kw)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"\n❌ 指令失敗：{' '.join(cmd)}\n{e.stderr}\n")
        sys.exit(1)


def parse_meme_file(meme_path: Path):
    txt = meme_path.read_text()
    first = txt.find("MOTIF ")
    header = txt[:first].strip().splitlines() if first != -1 else txt.strip().splitlines()
    motifs, evals = {}, {}
    for blk in txt.split("MOTIF ")[1:]:
        block = "MOTIF " + blk
        mid = blk.splitlines()[0].split()[0]
        m = re.search(r"\bE\s*=\s*([0-9.eE+-]+)", blk)
        e_val = float(m.group(1)) if m else float("inf")
        motifs[mid] = block
        evals[mid]  = e_val
    return header, motifs, evals


def filter_and_rename(raw_meme: Path, out_meme: Path,
                      e_thr: float, prefix: str):
    header, motifs, evals = parse_meme_file(raw_meme)
    with out_meme.open("w") as fh:
        fh.write("\n".join(header) + "\n\n")
        for mid, blk in motifs.items():
            if evals[mid] > e_thr:
                continue
            lines = blk.splitlines()
            parts = lines[0].split()
            parts[1] = f"{prefix}_{parts[1]}"
            lines[0] = " ".join(parts)
            fh.write("\n".join(lines).rstrip() + "\n\n")


def merge_meme(prepared_dir: Path, out_file: Path):
    files = sorted(prepared_dir.glob("*.meme"))
    if not files:
        sys.exit("❌ prepared 資料夾沒有 .meme 檔")
    head, *_ = parse_meme_file(files[0])
    out_file.write_text("\n".join(head) + "\n\n")
    with out_file.open("a") as fh:
        for f in files:
            _, motifs, _ = parse_meme_file(f)
            for blk in motifs.values():
                fh.write(blk.strip() + "\n\n")


def run_tomtom(meme_file: str, out_dir: str, thresh: float) -> str:
    run(["tomtom", "--oc", out_dir, "--thresh", str(thresh),
         "--dist", "pearson", meme_file, meme_file])
    return f"{out_dir}/tomtom.tsv"


def graph_dedupe(tsv: str, evals: dict):
    g = defaultdict(set)
    with open(tsv) as fh:
        for ln in fh:
            if (ln.startswith('#') or ln.startswith('Query_ID') or not ln.strip()):
                continue
            q, t, *_ = ln.split('\t')
            if q == t:
                continue
            g[q].add(t); g[t].add(q)

    visited, discard = set(), set()
    rep2dup = {}
    for node in g:
        if node in visited:
            continue
        q = deque([node]); comp = []
        while q:
            cur = q.popleft()
            if cur in visited:
                continue
            visited.add(cur)
            comp.append(cur)
            q.extend(g[cur] - visited)
        if len(comp) == 1:
            continue
        comp.sort(key=lambda m: (evals.get(m, float("inf")), m))
        rep, dups = comp[0], comp[1:]
        rep2dup[rep] = dups
        discard.update(dups)
    return rep2dup, discard


def write_outputs(header, motifs, evals, discard, rep2dup):
    kept = [m for m in motifs if m not in discard]
    kept.sort(key=lambda m: (evals.get(m, float("inf")), m))

    with open(FILTERED_MEME, "w") as fh:
        fh.write("\n".join(header) + "\n\n")
        for m in kept:
            fh.write(motifs[m].strip() + "\n\n")

    with open(KEPT_ID_TSV, "w") as fh:
        fh.write("experiment\tresult_dir\tmotif_id\te_value\n")
        for mid in kept:
            parts = mid.split("_")
            try:
                idx = parts.index("streme")
                exp = "_".join(parts[:idx])
                res = "_".join(parts[idx:idx+2])
                motif_id = "_".join(parts[idx+2:])
            except ValueError:
                exp, res, motif_id = "NA", "NA", mid
            fh.write(f"{exp}\t{res}\t{motif_id}\t{evals[mid]:.3g}\n")

    with open(REDUNDANT_TSV, "w") as fh:
        fh.write("representative\te_value\tdup_count\tduplicate_motifs\n")
        reps = sorted(rep2dup, key=lambda m: (evals.get(m, float('inf')), m))
        for r in reps:
            e_val = evals.get(r, float('nan'))
            fh.write(f"{r}\t{e_val:.3g}\t{len(rep2dup[r])}\t{';'.join(rep2dup[r])}\n")

    print(f"✅ Motif 去冗餘：總 {len(motifs)} → 保留 {len(kept)} (移除 {len(discard)})")
    print(f"➡ {FILTERED_MEME}\n➡ {KEPT_ID_TSV}\n➡ {REDUNDANT_TSV}")


# ╭──────────────── Motif-logo 共用函式 ─────────────────╮
def meme_info_matrix(pwm_df, nsites):
    entropy = -(pwm_df * np.log2(pwm_df.clip(lower=1e-9))).sum(axis=1)
    e_n = 3 / (2 * np.log(2) * nsites)
    return pwm_df.mul(2.0 - (entropy + e_n), axis=0).clip(lower=0)


def draw_logo_block(motif_block, title, save_path):
    rows = [list(map(float, ln.split()))
            for ln in motif_block.splitlines()
            if re.match(r'^\s*[0-9.]+\s', ln)]
    pwm = pd.DataFrame(rows, columns=['A','C','G','T'])
    m = re.search(r'nsites=\s*(\d+)', motif_block)
    nsites = int(m.group(1)) if m else 1
    info_mat = meme_info_matrix(pwm, nsites)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    lm.Logo(info_mat, ax=ax, color_scheme=COLOR_SCHEME)
    ax.set_title(title, fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    plt.tight_layout()
    fig.savefig(save_path, dpi=DPI, bbox_inches='tight')
    plt.close(fig)


# ╰───────────────────────────────────────────────────────╯
# ── logo-step ①：filtered.meme 前 N_LOGO_MOTIFS ──
def batch_plot_logos(filtered_meme: str,
                     out_dir: str,
                     n_motifs: int):
    txt = Path(filtered_meme).read_text()
    blocks = ["MOTIF "+b for b in txt.split("MOTIF ")[1:]]
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    for idx, blk in enumerate(blocks[:n_motifs], 1):
        motif_id = re.search(r'^MOTIF\s+(.+?)$', blk, re.MULTILINE).group(1).strip()
        safe = re.sub(r'[^A-Za-z0-9_-]', '_', motif_id)[:80]
        out_png = Path(out_dir) / f"{idx:02d}_{safe}.png"
        primary = motif_id.split()[0]
        try:
            species = primary.split('_')[1]
        except IndexError:
            species = "species?"
        motif_seq = primary.split('-')[-1]
        e_val = re.search(r'\bE\s*=\s*([0-9.eE+-]+)', blk)
        title = f"{species} motif: {motif_seq}"
        if e_val:
            title += f"  (E={float(e_val.group(1)):.2g})"
        draw_logo_block(blk, title, out_png)
        print(f"✔  {motif_seq} → {out_png}")
    print(f"✅ 前 {min(n_motifs,len(blocks))} 個 motif-logo 完成")


# ── logo-step ②：代表 motif（cluster reps）──
def parse_meme_blocks(meme_path):
    txt = Path(meme_path).read_text()
    mapping = {}
    for blk in txt.split("MOTIF ")[1:]:
        block = "MOTIF " + blk
        full_id = blk.splitlines()[0].strip()
        base_id = full_id.split(' STREME-')[0]
        mapping[base_id] = block
    return mapping


def batch_plot_rep_logos(filtered_meme: str,
                         redundant_tsv: str,
                         out_dir: str,
                         n_reps: int):
    df = pd.read_csv(redundant_tsv, sep='\t')
    df = df.dropna(subset=['representative'])
    df['e_value'] = pd.to_numeric(df['e_value'], errors='coerce')
    df.sort_values('e_value', inplace=True)
    reps = df.head(n_reps)

    id2block = parse_meme_blocks(filtered_meme)
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    for idx, row in enumerate(reps.itertuples(index=False), 1):
        rep_id  = row.representative
        dup_cnt = int(row.dup_count)
        blk = id2block.get(rep_id)
        if blk is None:
            print(f"⚠️  找不到 {rep_id} 的 MOTIF block，略過")
            continue

        safe = re.sub(r'[^A-Za-z0-9_-]', '_', rep_id)[:80]
        out_png = Path(out_dir) / f"{idx:02d}_{safe}.png"

        primary = rep_id.split('_')[-1]      # 最尾段含 motifSeq
        motif_seq = primary.split('-')[-1]
        species = rep_id.split('_')[1] if '_' in rep_id else 'species?'
        e_val = row.e_value
        title = f"{species} motif: {motif_seq} (E={e_val:.2g}, dup={dup_cnt})"
        draw_logo_block(blk, title, out_png)
        print(f"✔  {rep_id} → {out_png}")

    print(f"✅ 代表 motif-logo 完成，輸出於「{out_dir}/」")


# ────────────────────────── MAIN ──────────────────────────
def main():
    prepared_dir = Path(PREPARED_DIR)
    prepared_dir.mkdir(parents=True, exist_ok=True)
    Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

    # === 1. 轉檔 + 過濾 + 改名 ===
    n = 0
    for exp_dir in sorted(Path(DATA_DIR).iterdir()):
        if not exp_dir.is_dir():
            continue
        exp = exp_dir.name
        for streme_dir in sorted(d for d in exp_dir.iterdir()
                                 if d.is_dir() and d.name.startswith("streme_")):
            sub = streme_dir.name
            raw = None
            if (streme_dir / "streme.txt").is_file():
                raw = streme_dir / "streme.txt"
            elif (streme_dir / "streme.html").is_file():
                raw = streme_dir / "streme.html"
            else:
                for pat in ("*.meme", "*.txt", "*.html"):
                    files = list(streme_dir.glob(pat))
                    if files:
                        raw = files[0]
                        break
            if not raw:
                print(f"⚠️  找不到 STREME 結果：{streme_dir}")
                continue

            n += 1
            print(f"  [{n}] 處理 {exp}/{sub} → {raw.name}")

            tmp = prepared_dir / f"{exp}_{sub}.raw.meme"
            out = prepared_dir / f"{exp}_{sub}.meme"

            if raw.suffix.lower() == ".html":
                run(["meme2meme", str(raw)], stdout=tmp.open("w"))
            else:
                shutil.copy(raw, tmp)

            filter_and_rename(tmp, out, EVALUE_FILTER, f"{exp}_{sub}")
            tmp.unlink(missing_ok=True)

    if n == 0:
        sys.exit("❌ 未找到任何 STREME 資料夾，流程終止")

    # === 2. 合併 ===
    merge_meme(prepared_dir, Path(ALL_MEME))
    print("✔ 合併完成 → all.meme")

    # === 3. Tomtom 去冗餘 ===
    if Path(TEMP_DIR).exists():
        shutil.rmtree(TEMP_DIR, ignore_errors=True)

    header, motifs, evals = parse_meme_file(Path(ALL_MEME))
    if not motifs:
        print("⚠️  all.meme 無 motif，流程結束")
        return

    tsv = run_tomtom(ALL_MEME, TEMP_DIR, TOMTOM_THRESH)
    rep2dup, discard = graph_dedupe(tsv, evals)
    write_outputs(header, motifs, evals, discard, rep2dup)

    # === 4. filtered.meme 前 N_LOGO_MOTIFS ===
    batch_plot_logos(FILTERED_MEME, LOGO_OUT_DIR, N_LOGO_MOTIFS)

    # === 5. 代表 motif 前 N_REP_LOGOS ===
    if Path(REDUNDANT_TSV).is_file():
        batch_plot_rep_logos(FILTERED_MEME, REDUNDANT_TSV,
                             REP_LOGO_DIR, N_REP_LOGOS)
    else:
        print("⚠️  找不到 redundant_motif_reps.tsv，略過代表 motif-logo")


if __name__ == "__main__":
    main()
