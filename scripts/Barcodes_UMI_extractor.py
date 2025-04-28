import os, gzip, argparse, copy
import multiprocessing as mp
from functools import partial
from contextlib import contextmanager
from Bio import SeqIO
from Bio.Seq import Seq

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def levenshtein(a, b):
    if len(a) < len(b):
        return levenshtein(b, a)
    if len(b) == 0:
        return len(a)
    previous_row = range(len(b) + 1)
    for i, c1 in enumerate(a):
        current_row = [i + 1]
        for j, c2 in enumerate(b):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]

def get_header(seq, region):
    return seq[:region]

def get_tail(seq, region):
    return reverse_complement(seq[-region:])

def find_adapter(seq, adapter, max_dist=4):
    for i in range(len(seq) - len(adapter) + 1):
        window = seq[i:i+len(adapter)]
        dist = levenshtein(window.upper(), adapter.upper())
        if dist <= max_dist:
            return i + len(adapter)
    return -1

def has_polyA_or_polyT_tail(seq, threshold=8):
    return (
        any(seq[i:i+threshold].count("A") >= threshold for i in range(len(seq) - threshold + 1)) or
        any(seq[i:i+threshold].count("T") >= threshold for i in range(len(seq) - threshold + 1)))
    tail = seq[-scan_region:].upper()
    polyA_count = max(tail[i:i+threshold].count("A") for i in range(len(tail) - threshold + 1))
    polyT_count = max(tail[i:i+threshold].count("T") for i in range(len(tail) - threshold + 1))
    return polyA_count >= threshold or polyT_count >= threshold

@contextmanager
def poolcontext(*args, **kwargs):
    pool = mp.Pool(*args, **kwargs)
    yield pool
    pool.terminate()

def process_read(record, adapter, barcode_len, umi_len, min_len, scan_region, max_mismatch):
    record = copy.deepcopy(record)
    seq_str = str(record.seq)
    qual = record.letter_annotations.get("phred_quality", [])

    forward_seq = get_header(seq_str, scan_region)
    cut_pos = find_adapter(forward_seq, adapter, max_dist=max_mismatch)
    orientation = "forward"

    if cut_pos == -1:
        rev_seq = get_tail(seq_str, scan_region)
        cut_pos = find_adapter(rev_seq, adapter, max_dist=max_mismatch)
        orientation = "reverse"
        seq_str = reverse_complement(seq_str)
        qual = qual[::-1]

    if cut_pos != -1 and len(seq_str) >= cut_pos + barcode_len + umi_len:
        cb = seq_str[cut_pos:cut_pos + barcode_len]
        umi = seq_str[cut_pos + barcode_len:cut_pos + barcode_len + umi_len]
        cut_end = cut_pos + barcode_len + umi_len

        check_region = seq_str[cut_end:cut_end + 50]
        if has_polyA_or_polyT_tail(check_region):
            insert_seq = seq_str[cut_end:]
            insert_qual = qual[cut_end:]
            if len(insert_seq) == len(insert_qual) and len(insert_seq) >= min_len:
                record.letter_annotations.clear()
                record.seq = Seq(insert_seq)
                record.letter_annotations["phred_quality"] = insert_qual
                clean_id = record.id.split()[0]
                record.id = clean_id
                record.description = ""
                return record, [clean_id, orientation, cut_pos, cb, umi, len(insert_seq), round(sum(insert_qual)/len(insert_qual), 2)]
    return None, None

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, help="Input FASTQ file")
    parser.add_argument("-o", default="processed.fastq.gz", help="Output FASTQ file")
    parser.add_argument("-b", default="barcode_list.tsv.gz", help="Barcode info TSV file")
    parser.add_argument("-a", default="ACACGACGCTCTTCCGATCT", help="Adapter sequence")
    parser.add_argument("--barcode_len", type=int, default=38, help="Length of the barcode")
    parser.add_argument("--umi_len", type=int, default=8, help="Length of UMI")
    parser.add_argument("--min_read_length", type=int, default=200, help="Minimum read length after trimming")
    parser.add_argument("--scan_region", type=int, default=150, help="Region to scan for adapter")
    parser.add_argument("--ncores", type=int, default=4, help="Number of parallel workers")
    parser.add_argument("--max_mismatch", type=int, default=4, help="Maximum allowed mismatches for adapter detection")
    args = parser.parse_args()

    reader = gzip.open(args.i, "rt") if args.i.endswith(".gz") else open(args.i, "r")
    writer = gzip.open(args.o, "wt") if args.o.endswith(".gz") else open(args.o, "w")
    bc_out = gzip.open(args.b, "wt") if args.b.endswith(".gz") else open(args.b, "w")
    bc_out.write("\t".join(["read_id", "orientation", "BC_start", "BC", "UMI", "Seq_end", "mean_BC_quality"]) + "\n")


    records = list(SeqIO.parse(reader, "fastq"))

    count = 0
    adapter_found = 0
    adapter_not_found = 0
    failed_length = 0
    failed_polyA = 0
    failed_barcode = 0

    with poolcontext(processes=args.ncores) as pool:
        func = partial(process_read, adapter=args.a, barcode_len=args.barcode_len, umi_len=args.umi_len,
                      min_len=args.min_read_length, scan_region=args.scan_region, max_mismatch=args.max_mismatch)
        results = pool.map(func, records)

    for rec, info in results:
        if info:
            adapter_found += 1
            SeqIO.write(rec, writer, "fastq")
            bc_out.write("\t".join(str(x) for x in info) + "\n")
            count += 1
        else:
            adapter_not_found += 1
            seq_len = len(str(rec.seq)) if rec else 0
            if seq_len < args.barcode_len + args.umi_len + args.min_read_length:
                failed_length += 1
            elif rec and not has_polyA_or_polyT_tail(str(rec.seq)):
                failed_polyA += 1
            else:
                failed_barcode += 1

    reader.close()
    writer.close()
    bc_out.close()

    print(f"Total reads processed: {len(records)}")
    print(f"Reads kept: {count}")
    print(f"Reads with adapter found: {adapter_found}")
    print(f"Reads without adapter or failed QC: {adapter_not_found}")
    print(f" - Failed minimum length: {failed_length}")
    print(f" - Failed polyA/polyT check: {failed_polyA}")
    print(f" - Failed after adapter (barcode/UMI/polyA checks): {failed_barcode}")
