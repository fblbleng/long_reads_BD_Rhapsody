#!/usr/bin/env python3
"""
Barcode Assigner (Jake MAX FINAL Edition):
- Auto-tune merge_frac
- BK-Tree barcode merging
- Percentile UMI filtering
- Dual logging
- UMI diagnostics PDF
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import combinations
from rapidfuzz.distance import Levenshtein
from matplotlib.backends.backend_pdf import PdfPages
import logging
import time


# ─────────────────────────────────────────────────────────────────────────────
# BK-TREE CLASS
# ─────────────────────────────────────────────────────────────────────────────

class BKTree:
    def __init__(self, dist_fn):
        self.root = None
        self.dist_fn = dist_fn

    class Node:
        def __init__(self, word):
            self.word = word
            self.children = {}

    def add(self, word):
        if self.root is None:
            self.root = self.Node(word)
            return
        node = self.root
        while True:
            d = self.dist_fn(word, node.word)
            if d in node.children:
                node = node.children[d]
            else:
                node.children[d] = self.Node(word)
                break

    def search(self, word, max_dist):
        results = []
        if self.root is None:
            return results
        candidates = [self.root]
        while candidates:
            node = candidates.pop()
            d = self.dist_fn(word, node.word)
            if d <= max_dist:
                results.append(node.word)
            for k in node.children:
                if d - max_dist <= k <= d + max_dist:
                    candidates.append(node.children[k])
        return results


# ─────────────────────────────────────────────────────────────────────────────
# DATA PROCESSING + VISUALIZATION
# ─────────────────────────────────────────────────────────────────────────────

def count_umis_from_chunks(file_path, chunksize=500_000):
    umi_map = defaultdict(set)
    for chunk in pd.read_csv(file_path, sep='\t', chunksize=chunksize, usecols=['BC', 'UMI'], dtype=str):
        for bc, umi in zip(chunk['BC'], chunk['UMI']):
            umi_map[bc].add(umi)
    return pd.DataFrame({'BC': list(umi_map), 'UMI': [len(v) for v in umi_map.values()]})


def plot_umi_distribution(df, percentile_val=None):
    fig = plt.figure(figsize=(8, 5))
    plt.hist(df['UMI'], bins=100, color='orchid', edgecolor='black', log=True)
    plt.xlabel("UMI Count")
    plt.ylabel("Number of Barcodes (log)")
    plt.title("UMI Count Distribution")
    if percentile_val:
        plt.axvline(percentile_val, color='red', linestyle='--', label=f"≥ {int(percentile_val)} UMI (percentile)")
        plt.legend()
    plt.tight_layout()
    return fig

def plot_rank_curve(umi_counts, percentile_cutoff):
    sorted_counts = np.sort(umi_counts)[::-1]
    ranks = np.arange(1, len(sorted_counts)+1)
    fig = plt.figure(figsize=(8,5))
    plt.plot(ranks, sorted_counts, label='UMI counts', color='blue')
    plt.axhline(percentile_cutoff, color='red', linestyle='--', label=f'≥ {int(percentile_cutoff)} UMI')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Barcode Rank (log)")
    plt.ylabel("UMI Count (log)")
    plt.title("UMI Rank Curve")
    plt.legend()
    plt.tight_layout()
    return fig

def plot_percentile_curve(umi_counts):
    percentiles = np.arange(90, 100.1, 0.1)
    thresholds = [np.percentile(umi_counts, p) for p in percentiles]
    fig = plt.figure(figsize=(8,5))
    plt.plot(percentiles, thresholds, marker='o', color='darkgreen')
    plt.xlabel("Percentile")
    plt.ylabel("UMI Threshold")
    plt.title("UMI Threshold by Percentile")
    plt.grid(True)
    plt.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# MERGING + LOGGING
# ─────────────────────────────────────────────────────────────────────────────

def auto_tune_merge_frac(barcodes, percentile=95, sample_size=1000):
    sample = barcodes[:sample_size]
    dists = [Levenshtein.distance(a, b) for a, b in combinations(sample, 2)]
    est_ld = np.percentile(dists, percentile)
    bc_len = len(sample[0])
    merge_frac = est_ld / bc_len
    return merge_frac, dists


def merge_barcodes_bktree(cbs, max_dist, log):
    tree = BKTree(Levenshtein.distance)
    reps = []
    merged = {}
    for i, bc in enumerate(cbs):
        matches = tree.search(bc, max_dist)
        if matches:
            rep = matches[0]
            merged[rep].append(bc)
        else:
            tree.add(bc)
            reps.append(bc)
            merged[bc] = [bc]
        if i % 10000 == 0 and i > 0:
            log.info(f"🧬 Merged {i}/{len(cbs)} barcodes")
    return reps, merged


def setup_logging(log_file=None):
    log = logging.getLogger("BarcodeAssigner")
    log.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s', datefmt='%H:%M:%S')
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    log.addHandler(console_handler)
    if log_file:
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setFormatter(formatter)
        log.addHandler(file_handler)
    return log


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    start_time = time.time()

    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--input_raw_reads')
    g.add_argument('--input_counts')
    p.add_argument('-o', '--count_out', default='CB_counting_ext.tsv.gz')
    p.add_argument('-m', '--merged_out', default='CB_merged_list.txt')
    p.add_argument('--min_umi_per_cb', type=int, default=1)
    p.add_argument('--umi_percentile', type=float)
    p.add_argument('--merge_frac', type=str, default='0.2', help='"auto" or float')
    p.add_argument('--plot', action='store_true')
    p.add_argument('--log_file', help='Optional log file path')
    args = p.parse_args()

    log = setup_logging(args.log_file)
    log.info("🧬 Barcode Assigner STARTED")

    # ─ Load ─
    t1 = time.time()
    if args.input_raw_reads:
        df = count_umis_from_chunks(args.input_raw_reads)
        log.info("Loaded raw reads")
    else:
        df = pd.read_csv(args.input_counts, sep='\t', compression='infer', dtype={'BC': str})
        df['UMI'] = pd.to_numeric(df['UMI'], errors='coerce').fillna(0).astype(int)
        log.info("Loaded UMI count table")
    log.info(f"Loaded {len(df)} barcodes in {time.time() - t1:.2f}s")

    # ─ Filter ─
    df = df[df['UMI'] >= args.min_umi_per_cb].copy()
    df = df.sort_values('UMI', ascending=False).reset_index(drop=True)
    umi_counts = df['UMI'].values

    # ─ UMI Filtering ─
    if args.umi_percentile is not None:
        perc_val = np.percentile(umi_counts, args.umi_percentile)
        df['selected_raw'] = df['UMI'] >= perc_val
        log.info(f"Selected {df['selected_raw'].sum()} barcodes (UMI ≥ {int(perc_val)}, ≥ {args.umi_percentile}th percentile)")
    else:
        perc_val = np.percentile(umi_counts, 99)
        df['selected_raw'] = True
        log.info("No percentile filtering used; selected all")

    top_bcs = df[df['selected_raw']]['BC'].tolist()

    # ─ Plots to PDF ─
    if args.plot:
        fig1 = plot_umi_distribution(df, percentile_val=perc_val)
        fig2 = plot_rank_curve(umi_counts, perc_val)
        fig3 = plot_percentile_curve(umi_counts)
        with PdfPages("umi_analysis_plots.pdf") as pdf:
            pdf.savefig(fig1)
            pdf.savefig(fig2)
            pdf.savefig(fig3)
        log.info("📄 Saved UMI analysis plots → umi_analysis_plots.pdf")

    # ─ Merge Frac ─
    if args.merge_frac == 'auto':
        merge_frac, pre_dists = auto_tune_merge_frac(top_bcs)
        log.info(f"Auto-tuned merge_frac = {merge_frac:.3f}")
        if args.plot:
            plt.figure()
            plt.hist(pre_dists, bins=range(0, max(pre_dists)+2), color='skyblue', edgecolor='black')
            plt.title("Pre-Merge Pairwise LD")
            plt.tight_layout()
            plt.savefig("ld_premerge.png")
    else:
        merge_frac = float(args.merge_frac)

    max_dist = int(len(top_bcs[0]) * merge_frac)
    log.info(f"Using merge_frac={merge_frac:.3f} → LD ≤ {max_dist}")

    # ─ Merge ─
    t2 = time.time()
    reps, merged = merge_barcodes_bktree(top_bcs, max_dist, log)
    log.info(f"Merged {len(top_bcs)} → {len(reps)} barcodes in {time.time() - t2:.2f}s")

    # ─ Post-Merge LD plot ─
    if args.plot:
        post_dists = []
        for rep, group in merged.items():
            for bc in group:
                if bc != rep:
                    post_dists.append(Levenshtein.distance(bc, rep))
        plt.figure()
        plt.hist(post_dists, bins=range(0, max(post_dists)+2), color='orange', edgecolor='black')
        plt.title("Post-Merge Distances")
        plt.tight_layout()
        plt.savefig("ld_postmerge.png")

    # ─ Output ─
    rep_map = {bc: rep for rep, members in merged.items() for bc in members}
    df['rep_id'] = df['BC'].map(rep_map).fillna(df['BC'])

    df.to_csv(args.count_out, sep='\t', index=False, compression='gzip',
              columns=['BC', 'UMI', 'rep_id', 'selected_raw'])
    log.info(f"Wrote count table → {args.count_out}")

    with open(args.merged_out, 'w') as fw:
        for rep in reps:
            fw.write(f"{rep}:\t{','.join(merged[rep])}\n")
    log.info(f"Wrote merged list → {args.merged_out}")

    log.info(f"✅ Pipeline complete in {time.time() - start_time:.2f}s")


if __name__ == '__main__':
    main()
