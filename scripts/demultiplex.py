#!/usr/bin/env python3
"""
Demultiplexer with % progress and fuzzy matching (BK-tree)
"""

import sys
import gzip
import argparse
import itertools
from Bio import SeqIO
from multiprocessing import Pool
from Levenshtein import distance as ld
from pybktree import BKTree
from functools import lru_cache
import time
import threading

# Globals for workers
READ_MAP = {}
REP_LIST = []
MAX_MM = 0
REP_MAP_DICT = {}
BK_TREE = None

progress_lock = threading.Lock()
reads_done = 0
total_reads = 1


def load_barcode_table(fn):
    mapping = {}
    opener = gzip.open if fn.endswith('.gz') else open
    with opener(fn, 'rt') as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if cols[0].lower() in ('read_id', 'readid'):
                continue
            rid, cb, umi = cols[0], cols[3], cols[4]
            mapping[rid] = (cb, umi)
    return mapping


def load_merged_list(fn):
    rep_map = {}
    with open(fn) as fh:
        for line in fh:
            if not line:
                continue
            rep, members = line.strip().split(':', 1)
            for m in members.split(','):
                rep_map[m] = rep
    return rep_map


@lru_cache(maxsize=100000)
def fuzzy_match_cb(raw_cb):
    matches = BK_TREE.find(raw_cb, MAX_MM)
    if matches:
        return sorted(matches, key=lambda x: x[1])[0][0]
    return raw_cb


def init_worker(rmap, rep_list, max_mm, rep_map_dict):
    global READ_MAP, REP_LIST, MAX_MM, REP_MAP_DICT, BK_TREE
    READ_MAP = rmap
    REP_LIST = rep_list
    MAX_MM = max_mm
    REP_MAP_DICT = rep_map_dict
    if MAX_MM > 0:
        BK_TREE = BKTree(ld, REP_LIST)


def annotate_batch(batch):
    out = []
    for rec in batch:
        rid = rec.id.split()[0]
        info = READ_MAP.get(rid)
        if not info:
            continue
        raw_cb, umi = info
        rep = REP_MAP_DICT.get(raw_cb)
        if rep is None and MAX_MM > 0:
            rep = fuzzy_match_cb(raw_cb)
        elif rep is None:
            rep = raw_cb
        rec.id = f"{rid}_CB:{raw_cb}_UMI:{umi}_REP:{rep}"
        rec.description = ''
        out.append(rec)
    return out


def chunkify(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def count_reads_in_file(fastq_path):
    opener = gzip.open if fastq_path.endswith('.gz') else open
    with opener(fastq_path, 'rt') as f:
        return sum(1 for _ in f) // 4


def main():
    parser = argparse.ArgumentParser(description="Demultiplexer with % progress")
    parser.add_argument('-i', '--input_fastq', required=True, help='Input FASTQ (file only, not -)')
    parser.add_argument('-t', '--table', required=True, help='barcode_list.tsv(.gz)')
    parser.add_argument('-m', '--merged_list', required=True, help='CB_merged_list.txt')
    parser.add_argument('-o', '--output_fastq', default='-', help='Output FASTQ (- for stdout)')
    parser.add_argument('--max_mismatch', type=int, default=0, help='Max Levenshtein distance')
    parser.add_argument('-T', '--threads', type=int, default=4, help='Worker count')
    parser.add_argument('--chunksize', type=int, default=10000, help='Reads per batch')
    args = parser.parse_args()

    if args.input_fastq == "-":
        sys.exit("❌ Percentage mode requires a real input file, not stdin (-i -)")

    # Pre-count reads
    global total_reads
    print(f"🧮 Counting total reads in {args.input_fastq}...", file=sys.stderr)
    total_reads = count_reads_in_file(args.input_fastq)
    print(f"✅ Total reads: {total_reads:,}", file=sys.stderr)

    read_map = load_barcode_table(args.table)
    rep_map = load_merged_list(args.merged_list)
    rep_list = list(set(rep_map.values()))

    fin = gzip.open(args.input_fastq, 'rt') if args.input_fastq.endswith('.gz') else open(args.input_fastq)
    fout = gzip.open(args.output_fastq, 'wt') if args.output_fastq != '-' else sys.stdout

    pool = Pool(processes=args.threads,
                initializer=init_worker,
                initargs=(read_map, rep_list, args.max_mismatch, rep_map))

    reader = SeqIO.parse(fin, 'fastq')
    start_time = time.time()

    for batch in iter(lambda: list(itertools.islice(reader, args.chunksize)), []):
        sub_batches = list(chunkify(batch, max(1, len(batch) // args.threads)))
        results = pool.map(annotate_batch, sub_batches)

        for recs in results:
            SeqIO.write(recs, fout, 'fastq')

        with progress_lock:
            global reads_done
            reads_done += len(batch)
            elapsed = time.time() - start_time
            percent = (reads_done / total_reads) * 100
            rate = reads_done / elapsed if elapsed > 0 else 0
            print(f"[{time.strftime('%H:%M:%S')}] {reads_done:,}/{total_reads:,} "
                  f"({percent:.2f}%) — {rate:,.0f} reads/sec", file=sys.stderr)

    pool.close()
    pool.join()


if __name__ == '__main__':
    main()
