#!/usr/bin/env python3
"""
Full Assigner Pipeline: Knee detection, dynamic rescue, and barcode merging
"""
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- Levenshtein distance ---
def levenshtein(a, b):
    if len(a) < len(b):
        return levenshtein(b, a)
    if not b:
        return len(a)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        curr = [i]
        for j, cb in enumerate(b, 1):
            cost = 0 if ca == cb else 1
            curr.append(min(curr[j-1] + 1,
                            prev[j]   + 1,
                            prev[j-1] + cost))
        prev = curr
    return prev[-1]

# --- Merge barcodes within edit distance ---
def merge_barcodes(cbs, max_dist):
    reps = []
    merged = {}
    for bc in cbs:
        placed = False
        for rep in reps:
            if levenshtein(bc, rep) <= max_dist:
                merged[rep].append(bc)
                placed = True
                break
        if not placed:
            reps.append(bc)
            merged[bc] = [bc]
    return reps, merged

# --- Knee & Rescue detection ---
def detect_knee_and_rescue(umi_counts, smooth_res, rescue_frac):
    n = len(umi_counts)
    if n < 2:
        return 0, 0
    ranks = np.arange(1, n+1)
    log_counts = np.log10(umi_counts + 1)
    slopes = np.gradient(log_counts)
    log_ranks = np.log10(ranks)
    bin_idx = np.floor(log_ranks / smooth_res).astype(int)
    df = pd.DataFrame({'bin': bin_idx, 'slope': slopes})
    medians = df.groupby('bin')['slope'].median()
    best_bin = medians.idxmin()
    idxs = df.index[df['bin']==best_bin].tolist()
    if not idxs:
        knee = int(np.argmin(slopes))
    else:
        knee = idxs[0]
    rescue = min(int(knee * (1 + rescue_frac)), n-1)
    return knee, rescue

def main():
    p = argparse.ArgumentParser(
        description="Assigner: knee detection + dynamic rescue + barcode merging"
    )
    p.add_argument('-i','--input', required=True,
                   help='CB count table (BC and UMI columns, tab-delimited)')
    p.add_argument('-o','--count_out', default='CB_counting_ext.tsv.gz',
                   help='Extended UMI count table output')
    p.add_argument('-m','--merged_out', default='CB_merged_list.txt',
                   help='Merged barcode list output')
    p.add_argument('--min_umi_per_cb', type=int, default=1,
                   help='Filter barcodes with fewer than this UMI count')
    p.add_argument('--smooth_res', type=float, default=0.001,
                   help='Log10‐rank bin width for slope smoothing')
    p.add_argument('--rescue_frac', type=float, default=0.10,
                   help='Fractional extension past knee (e.g. 0.1)')
    p.add_argument('--merge_dist', type=int, default=None,
                   help='(Optional) absolute LD threshold to merge barcodes')
    p.add_argument('--merge_frac', type=float, default=None,
                   help='(Optional) fraction of BC length to use as LD merge threshold (e.g. 0.2)')
    p.add_argument('--forced_no', type=int, default=0,
                   help='If >0, force fixed number of barcodes (skip knee)')
    p.add_argument('--plot', action='store_true',
                   help='Save knee+rescue plot to "knee_rescue.png"')
    args = p.parse_args()

    # --- 1) Load & filter ---
    df = pd.read_csv(args.input, sep='\t', compression='infer',
                     dtype={'BC': str})
    df['UMI'] = pd.to_numeric(df['UMI'], errors='coerce').fillna(0).astype(int)
    df = df[df['UMI'] >= args.min_umi_per_cb].copy()
    df = df.sort_values('UMI', ascending=False).reset_index(drop=True)

    umi_counts = df['UMI'].values

    # --- 2) Knee detection ---
    if args.forced_no > 0:
        knee = rescue = args.forced_no - 1
        print(f"Forced to top {args.forced_no} barcodes → knee={knee+1}")
    else:
        knee, rescue = detect_knee_and_rescue(
            umi_counts, args.smooth_res, args.rescue_frac
        )
        print(f"Knee at rank {knee+1}, rescue at {rescue+1}")

    # mark selected
    df['selected_raw'] = False
    df.loc[:rescue, 'selected_raw'] = True

    # --- 3) Determine merge threshold ---
    top_bcs = df.loc[df['selected_raw'], 'BC'].tolist()
    if args.merge_dist is not None:
        max_dist = args.merge_dist
    elif args.merge_frac is not None:
        # use length of first barcode as reference
        bc_len = len(top_bcs[0]) if top_bcs else 0
        max_dist = int(bc_len * args.merge_frac)
    else:
        max_dist = 1

    # --- 4) Merge barcodes ---
    reps, merged = merge_barcodes(top_bcs, max_dist)
    print(f"Collapsed {len(top_bcs)} → {len(reps)} representative barcodes (LD≤{max_dist})")

    # annotate rep IDs in df
    rep_map = {}
    for rep, members in merged.items():
        for bc in members:
            rep_map[bc] = rep
    df['rep_id'] = df['BC'].map(rep_map).fillna(df['BC'])

    # --- 5) Write extended count table ---
    out_cols = ['BC','UMI','rep_id','selected_raw']
    df.to_csv(args.count_out, sep='\t', index=False, compression='gzip',
              columns=out_cols)

    # --- 6) Write merged list ---
    with open(args.merged_out, 'w') as fw:
        for rep in reps:
            fw.write(f"{rep}:\t{','.join(merged[rep])}\n")

    # --- 7) Optional plot ---
    if args.plot:
        import matplotlib.pyplot as plt
        ranks = np.arange(len(umi_counts))
        log_counts = np.log10(umi_counts + 1)
        plt.figure()
        plt.plot(ranks, log_counts, label='log10(UMI+1)')
        plt.axvline(knee, color='r', linestyle='--', label='Knee')
        plt.axvline(rescue, color='g', linestyle='--', label='Rescue')
        plt.xlabel('Rank')
        plt.ylabel('log10(UMI+1)')
        plt.title('Knee & Rescue')
        plt.legend()
        plt.savefig('knee_rescue.png')
        print("Saved plot to knee_rescue.png")

if __name__ == '__main__':
    main()

