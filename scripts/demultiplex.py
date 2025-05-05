#!/usr/bin/env python3
import os
import gzip
import argparse
from functools import partial
from multiprocessing import Pool
from Bio import SeqIO
import pandas as pd

# globals filled in initializer
READ_MAP = None      # read_id -> {'BC':…, 'UMI':…}
REPS     = None      # list of representative barcodes
RAW2REP  = None      # raw_BC -> representative_CB
MAX_MM   = None      # max mismatch

def levenshtein(a: str, b: str) -> int:
    if len(a) < len(b):
        return levenshtein(b, a)
    if not b:
        return len(a)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        curr = [i]
        for j, cb in enumerate(b, 1):
            cost = 0 if ca == cb else 1
            curr.append(min(curr[j-1] + 1, prev[j] + 1, prev[j-1] + cost))
        prev = curr
    return prev[-1]

def init_worker(read_map, reps, raw2rep, max_mm):
    global READ_MAP, REPS, RAW2REP, MAX_MM
    READ_MAP = read_map
    REPS     = reps
    RAW2REP  = raw2rep
    MAX_MM   = max_mm

def worker(rec):
    """
    Annotate one SeqRecord: return (rep, fastq_str) or None if drop.
    """
    rid = rec.id.split('_CB:')[0]
    info = READ_MAP.get(rid)
    if not info:
        return None
    raw_bc, umi = info['BC'], info['UMI']
    rep = RAW2REP.get(raw_bc)
    if rep is None and MAX_MM > 0:
        # fuzzy‐match
        best, bd = None, MAX_MM + 1
        for cand in REPS:
            d = levenshtein(raw_bc, cand)
            if d < bd:
                bd, best = d, cand
                if d == 0:
                    break
        if bd <= MAX_MM:
            rep = best
    if not rep:
        return None

    # build new header
    new_id = f"{rid}_CB:{raw_bc}_UMI:{umi}_REP:{rep}"
    rec.id = new_id
    rec.description = ''
    fastq_str = rec.format('fastq')
    return rep, fastq_str

def demultiplex(fastq_in, output_dir, threads, chunk_size=5000):
    os.makedirs(output_dir, exist_ok=True)
    writers = {}
    total = kept = 0

    opener = gzip.open if fastq_in.endswith('.gz') else open
    with opener(fastq_in, 'rt') as fin, Pool(threads, initializer=init_worker,
                                             initargs=(READ_MAP, REPS, RAW2REP, MAX_MM)) as pool:
        it = SeqIO.parse(fin, 'fastq')
        while True:
            chunk = list()
            try:
                for _ in range(chunk_size):
                    chunk.append(next(it))
            except StopIteration:
                pass
            if not chunk:
                break

            # annotate in parallel
            for res in pool.map(worker, chunk):
                total += 1
                if res:
                    rep, fq = res
                    if rep not in writers:
                        path = os.path.join(output_dir, f"{rep}.fastq.gz")
                        writers[rep] = gzip.open(path, 'wt')
                    writers[rep].write(fq)
                    kept += 1

    # close all
    for fh in writers.values():
        fh.close()

    print(f"Total reads scanned: {total}")
    print(f"Total reads kept:  {kept}")
    print(f"Wrote {len(writers)} per-cell FASTQ files in '{output_dir}'")

def load_read_mapping(map_file):
    df = pd.read_csv(map_file, sep='\t', compression='infer',
                     usecols=['read_id','BC','UMI'], dtype=str)
    df = df.dropna(subset=['read_id','BC','UMI'])
    return df.set_index('read_id')[['BC','UMI']].to_dict('index')

def load_merged_list(rep_file):
    raw2rep = {}
    reps = []
    with open(rep_file) as f:
        for line in f:
            rep, members = line.strip().split(':',1)
            rep = rep.strip()
            reps.append(rep)
            for bc in members.strip().split(','):
                raw2rep[bc] = rep
    return sorted(set(reps)), raw2rep

if __name__ == '__main__':
    p = argparse.ArgumentParser(
        description="Parallel demultiplex trimmed FASTQ by Representative_CB"
    )
    p.add_argument('-i','--input_fastq', required=True,
                   help='Trimmed FASTQ (processed.fastq.gz)')
    p.add_argument('-m','--mapping', required=True,
                   help='barcode_list.tsv.gz (read_id,BC,UMI)')
    p.add_argument('-r','--rep_file', required=True,
                   help='CB_merged_list.txt (rep:member,member,...)')
    p.add_argument('-o','--output_dir', default='demux_fastq',
                   help='Output directory for per-cell FASTQ')
    p.add_argument('--max_mismatch', type=int, default=1,
                   help='Max LD to fuzzy‐map raw BC→rep')
    p.add_argument('-t','--threads', type=int, default=4,
                   help='Number of worker processes')
    args = p.parse_args()

    print("Loading read→BC,UMI mapping…")
    read_map = load_read_mapping(args.mapping)
    print(f" ↳ {len(read_map)} reads mapped")

    print("Loading rep CB list…")
    reps, raw2rep = load_merged_list(args.rep_file)
    print(f" ↳ {len(reps)} representative barcodes")

    # set globals for initializer
    READ_MAP = read_map
    REPS     = reps
    RAW2REP  = raw2rep
    MAX_MM   = args.max_mismatch

    print("Demultiplexing in parallel…")
    demultiplex(args.input_fastq, args.output_dir, args.threads)
    print("All done.")
