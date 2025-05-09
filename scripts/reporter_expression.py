#!/usr/bin/env python3
"""
reporter_expression.py

Quantify gene expression per cell from curated BAMs using FeatureCounts.
"""

import os
import sys
import argparse
import subprocess
import pandas as pd

def run(cmd):
    """Run shell command, exit on error."""
    print(f"> {cmd}")
    proc = subprocess.Popen(cmd, shell=True)
    proc.communicate()
    if proc.returncode != 0:
        sys.stderr.write(f"❌ Command failed (exit {proc.returncode}):\n  {cmd}\n")
        sys.exit(proc.returncode)

def main():
    p = argparse.ArgumentParser(
        description="Generate gene×cell UMI count matrix from per-cell curated BAMs"
    )
    p.add_argument('-d','--workdir', required=True,
                   help='Directory containing *.curated.bam (+ .bai) for each cell')
    p.add_argument('--gtf', required=True,
                   help='Gene annotation GTF file')
    p.add_argument('-o','--out', default='matrix.tsv',
                   help='Output counts matrix (TSV)')
    p.add_argument('--sel_bc_o', default=None,
                   help='Optional: file of selected barcodes (one per line) to subset cells')
    p.add_argument('-t','--threads', type=int, default=4,
                   help='Number of threads for FeatureCounts')
    p.add_argument('--featurecounts', default='featureCounts',
                   help='Path to featureCounts binary')
    p.add_argument('--keep_tmp', action='store_true',
                   help='If set, keep intermediate featureCounts output')
    args = p.parse_args()

    # 1) discover all curated BAMs
    bam_files = []
    for fn in sorted(os.listdir(args.workdir)):
        if fn.endswith('.curated.bam'):
            cell = fn[:-len('.curated.bam')]
            bam_files.append((cell, os.path.join(args.workdir, fn)))

    if not bam_files:
        sys.exit("❌ No *.curated.bam files found in " + args.workdir)

    # 2) optionally filter by sel_bc_o
    if args.sel_bc_o:
        sel = set(open(args.sel_bc_o).read().split())
        bam_files = [(c,p) for c,p in bam_files if c in sel]
        if not bam_files:
            sys.exit("❌ After subsetting, no cells remain.")

    cells, paths = zip(*bam_files)

    # 3) run featureCounts
    tmp_counts = "featureCounts.raw.txt"
    cmd = (
        f"{args.featurecounts} "
        f"-T {args.threads} "
        f"-a {args.gtf} "
        f"-o {tmp_counts} "
        + " ".join(paths)
    )
    run(cmd)

    # 4) parse featureCounts output
    # skip first line (header), read until you hit columns
    df = pd.read_csv(tmp_counts, sep='\t', comment='#', header=0)
    # columns: Geneid, Chr, Start, End, Strand, Length, <cell1.bam>, <cell2.bam>, ...
    count_cols = [c for c in df.columns if c.endswith('.bam')]
    # strip the .curated.bam suffix to get cell IDs
    new_cols = {col: os.path.basename(col).replace('.curated.bam','') for col in count_cols}
    df = df[['Geneid'] + count_cols].rename(columns=new_cols)

    # 5) write matrix.tsv
    df.to_csv(args.out, sep='\t', index=False)
    print(f"✅ Wrote gene×cell matrix to {args.out}")

    # 6) cleanup
    if not args.keep_tmp:
        os.remove(tmp_counts)
        for ext in ('.counts', '.featureCounts', '.summary'):
            fn = tmp_counts + ext
            if os.path.exists(fn):
                os.remove(fn)

if __name__ == '__main__':
    main()
