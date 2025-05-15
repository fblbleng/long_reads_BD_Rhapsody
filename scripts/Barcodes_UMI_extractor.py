#!/usr/bin/env python3
import os, gzip, argparse, itertools, multiprocessing as mp
from functools import partial
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm           #  ← NEW

# ---------- helper functions (unchanged) ----------
def reverse_complement(s): return str(Seq(s).reverse_complement())
def levenshtein(a, b):
    if len(a) < len(b): a, b = b, a
    prev = range(len(b) + 1)
    for i, c1 in enumerate(a):
        curr = [i + 1]
        for j, c2 in enumerate(b):
            curr.append(min(prev[j+1]+1, curr[j]+1, prev[j] + (c1 != c2)))
        prev = curr
    return prev[-1]

def get_header(s, n): return s[:n]
def get_tail(s, n):  return reverse_complement(s[-n:])

def find_adapter(seq, adapter, max_dist=4):
    L = len(adapter)
    for i in range(len(seq) - L + 1):
        if levenshtein(seq[i:i+L].upper(), adapter.upper()) <= max_dist:
            return i + L
    return -1

def has_polyA_or_polyT_tail(seq, threshold=8):
    for i in range(len(seq) - threshold + 1):
        w = seq[i:i+threshold]
        if w.count("A") >= threshold or w.count("T") >= threshold:
            return True
    return False

def process_read(record, adapter, barcode_len, umi_len,
                 min_len, scan_region, max_mismatch):
    seq = str(record.seq)
    qual = record.letter_annotations["phred_quality"]

    cut = find_adapter(get_header(seq, scan_region), adapter, max_mismatch)
    orientation = "forward"
    if cut == -1:
        cut = find_adapter(get_tail(seq, scan_region), adapter, max_mismatch)
        if cut == -1:
            return None
        seq = reverse_complement(seq)
        qual = qual[::-1]
        orientation = "reverse"

    if len(seq) < cut + barcode_len + umi_len:
        return None

    cb  = seq[cut:cut+barcode_len]
    umi = seq[cut+barcode_len : cut+barcode_len+umi_len]
    insert_seq  = seq[cut+barcode_len+umi_len:]
    insert_qual = qual[cut+barcode_len+umi_len:]

    if len(insert_seq) < min_len or not has_polyA_or_polyT_tail(insert_seq[:50]):
        return None

    rec_out = SeqIO.SeqRecord(Seq(insert_seq),
                              id=record.id.split()[0],
                              description="")
    rec_out.letter_annotations["phred_quality"] = insert_qual
    info = [rec_out.id, orientation, cut, cb, umi,
            len(insert_seq), round(sum(insert_qual)/len(insert_qual), 2)]
    return rec_out, info

# ---------- streaming driver with progress bar ----------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, help="Input FASTQ(.gz)")
    parser.add_argument("-o", default="processed.fastq.gz")
    parser.add_argument("-b", default="barcode_list.tsv.gz")
    parser.add_argument("-a", default="ACACGACGCTCTTCCGATCT")
    parser.add_argument("--barcode_len", type=int, default=38)
    parser.add_argument("--umi_len",     type=int, default=8)
    parser.add_argument("--min_read_length", type=int, default=200)
    parser.add_argument("--scan_region",     type=int, default=150)
    parser.add_argument("--ncores", type=int, default=4)
    parser.add_argument("--chunksize", type=int, default=10000)
    parser.add_argument("--max_mismatch", type=int, default=4)
    args = parser.parse_args()

    reader = gzip.open(args.i, "rt") if args.i.endswith(".gz") else open(args.i)
    writer = gzip.open(args.o, "wt") if args.o.endswith(".gz") else open(args.o, "w")
    bc_out = gzip.open(args.b, "wt") if args.b.endswith(".gz") else open(args.b, "w")
    bc_out.write("\t".join(["read_id","orientation","BC_start","BC","UMI",
                            "Seq_end","mean_BC_quality"]) + "\n")

    record_iter = SeqIO.parse(reader, "fastq")
    func = partial(process_read, adapter=args.a, barcode_len=args.barcode_len,
                   umi_len=args.umi_len, min_len=args.min_read_length,
                   scan_region=args.scan_region, max_mismatch=args.max_mismatch)

    kept = total = 0
    with mp.Pool(args.ncores) as pool, tqdm(total=None, desc="reads", unit="read") as pbar:
        batched = iter(lambda: list(itertools.islice(record_iter, args.chunksize)), [])
        for batch in batched:
            for rec_info in pool.imap_unordered(func, batch, chunksize=100):
                total += 1
                pbar.update()          #  ← NEW
                if rec_info:
                    rec, info = rec_info
                    SeqIO.write(rec, writer, "fastq")
                    bc_out.write("\t".join(map(str, info)) + "\n")
                    kept += 1

    reader.close(); writer.close(); bc_out.close()
    print(f"\nProcessed {total:,} reads  •  kept {kept:,} ({kept/total:.2%})")

if __name__ == "__main__":
    main()
