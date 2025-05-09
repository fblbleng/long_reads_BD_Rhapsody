#!/usr/bin/env python3
"""
curator.py

Curator Pipeline: per-cell mapping, soft-clip filtering, UMI collapsing (SPOA), remapping & final BAM.
"""
import os
import sys
import time
import gzip
import argparse
import subprocess
import multiprocessing as mp
from functools import partial
from contextlib import contextmanager

import pysam
from Bio import SeqIO

# -----------------------------------------------------------------------------
def sys_run(cmd):
    """Run shell command, return (returncode, stdout, stderr)."""
    proc = subprocess.Popen(cmd, shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    out, err = proc.communicate()
    return proc.returncode, out, err

def check_mapping_db(args):
    """Return the reference argument for minimap2."""
    return args.idx_genome if args.idx_genome else args.ref_genome

# -----------------------------------------------------------------------------
def filter_softclipping(thr, in_bam, out_bam, drop_bam, samtools):
    """
    Keep reads with (matches)/(matches+softclips) ≥ thr;
    drop unmapped or fully‐softclipped reads.
    """
    bam_i = pysam.AlignmentFile(in_bam, "rb")
    bam_o = pysam.AlignmentFile(out_bam, "wb", template=bam_i)
    dropped = set()

    for r in bam_i.fetch(until_eof=True):
        if r.is_unmapped or r.cigartuples is None:
            dropped.add(r.query_name)
            continue
        m = sum(length for op, length in r.cigartuples if op == 0)
        s = sum(length for op, length in r.cigartuples if op in (4, 5))
        if m + s > 0 and (m / (m + s)) >= thr:
            bam_o.write(r)
        else:
            dropped.add(r.query_name)

    bam_i.close()
    bam_o.close()

    bam_i = pysam.AlignmentFile(in_bam, "rb")
    bam_d = pysam.AlignmentFile(drop_bam, "wb", template=bam_i)
    for r in bam_i.fetch(until_eof=True):
        if r.query_name in dropped:
            bam_d.write(r)
    bam_i.close()
    bam_d.close()

# -----------------------------------------------------------------------------
def split_reads_by_pos(bam, coord_buffer):
    """
    Group reads into clusters by chromosome and start position drift ≤ coord_buffer.
    Each cluster is a list of dicts: {'query_name','chr_name','start'}.
    """
    clusters, buf = [], []
    prev_chr, prev_pos = None, None
    for r in bam.fetch(until_eof=True):
        if prev_chr is None:
            prev_chr, prev_pos = r.reference_name, r.reference_start
        if r.reference_name != prev_chr or (r.reference_start - prev_pos) > coord_buffer:
            clusters.append(buf)
            buf = []
        buf.append({
            'query_name': r.query_name,
            'chr_name':   r.reference_name,
            'start':      r.reference_start
        })
        prev_chr, prev_pos = r.reference_name, r.reference_start
    if buf:
        clusters.append(buf)
    return clusters

# -----------------------------------------------------------------------------
def build_read_seq_dict(fq_path):
    """
    Build dict: read_id → sequence, from FASTQ (handles .gz).
    """
    d = {}
    opener = gzip.open if fq_path.endswith('.gz') else open
    with opener(fq_path, 'rt') as fin:
        for rec in SeqIO.parse(fin, 'fastq'):
            d[rec.id] = str(rec.seq)
    return d

# -----------------------------------------------------------------------------
def levenshtein(a, b):
    """Compute Levenshtein distance."""
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

# -----------------------------------------------------------------------------
def collapse_umi(cluster, seq_dict,
                 umi_ld, coord_buffer,
                 max_dup, spoa_path, tmp_dir):
    """
    Collapse UMIs within a genomic cluster.
    - cluster: list of {'query_name','chr_name','start'}
    - seq_dict: read_id → sequence map
    - umi_ld: Levenshtein distance threshold for merging UMIs
    - coord_buffer: allowed positional drift (unused here, but kept for signature)
    - max_dup: max number of reads per UMI to feed to SPOA
    - spoa_path: full path to SPOA binary
    - tmp_dir: where to write temp FASTA
    Returns (consensus_dict, duplicates_set, singletons_list).
    """
    # 1) group by UMI tag in read name
    umi_groups = {}
    for r in cluster:
        rid = r['query_name']
        if '_UMI:' not in rid:
            continue
        umi = rid.rsplit('_UMI:',1)[-1]
        umi_groups.setdefault(umi, []).append(r)

    # 2) merge UMIs by LD ≤ umi_ld
    merged = {}
    seen = set()
    for umi in umi_groups:
        if umi in seen:
            continue
        group = [umi]
        seen.add(umi)
        for umi2 in umi_groups:
            if umi2 not in seen and levenshtein(umi, umi2) <= umi_ld:
                group.append(umi2)
                seen.add(umi2)
        merged[umi] = group

    consensus  = {}
    duplicates  = set()
    singletons  = []

    # 3) for each merged UMI group, build consensus if >1
    for rep_umi, members in merged.items():
        reads = []
        for u in members:
            reads.extend(umi_groups[u])
        if len(reads) > 1:
            # write temp FASTA
            fa_path = os.path.join(tmp_dir, f"{rep_umi}.fa")
            with open(fa_path, 'wt') as fw:
                for idx, r in enumerate(reads):
                    if idx >= max_dup:
                        break
                    rid = r['query_name']
                    seq = seq_dict.get(rid, "")
                    fw.write(f">{rid}\n{seq}\n")
                    duplicates.add(rid)
            # run SPOA
            rc, out, err = sys_run(f"{spoa_path} {fa_path}")
            os.remove(fa_path)
            if rc == 0:
                lines = out.decode().splitlines()
                if len(lines) >= 2:
                    consensus[reads[0]['query_name']] = lines[1]
        else:
            # singleton UMI: keep the single read
            singletons.append(reads[0]['query_name'])

    return consensus, duplicates, singletons

# -----------------------------------------------------------------------------
def curate_one(fq_path, options):
    """Process one cell FASTQ end-to-end."""
    cell = os.path.basename(fq_path).split('.')[0]
    tmp_prefix = os.path.join(options.tmp_dir, cell)
    os.makedirs(options.tmp_dir, exist_ok=True)
    db = check_mapping_db(options)

    # 1) map → unsorted.bam
    uns = f"{tmp_prefix}.unsorted.bam"
    cmd = (f"{options.minimap2} -ax splice {db} {fq_path} "
           f"| {options.samtools} view -Sb - > {uns}")
    rc, _, _ = sys_run(cmd)
    if rc:
        sys.stderr.write(f"❌ minimap2 failed on {cell}\n")
        return None

    # 2) filter soft clips
    filt = f"{tmp_prefix}.filtered.bam"
    drop = f"{tmp_prefix}.dropped.bam"
    filter_softclipping(options.softclipping_thr, uns, filt, drop, options.samtools)

    # 3) index filtered BAM
    sys_run(f"{options.samtools} index {filt}")

    # 4) load original sequences
    seq_dict = build_read_seq_dict(fq_path)

    # 5) split reads by position clusters
    bam_i = pysam.AlignmentFile(filt, "rb")
    clusters = split_reads_by_pos(bam_i, options.coord_buffer)
    bam_i.close()

    # 6) collapse UMIs in each cluster
    all_cons, all_dups, all_sins = {}, set(), []
    for cl in clusters:
        c, d, s = collapse_umi(
            cl, seq_dict,
            options.umi_ld,
            options.coord_buffer,
            options.max_umi_duplicates,
            options.spoa,
            options.tmp_dir
        )
        all_cons.update(c)
        all_dups |= d
        all_sins += s

    # 7) write consensus FASTA
    fa_cons = f"{tmp_prefix}.consensus.fa"
    with open(fa_cons, 'wt') as fw:
        for rid, seq in all_cons.items():
            fw.write(f">{rid}\n{seq}\n")

    # 8) write singleton BAM
    sin = f"{tmp_prefix}.singleton.bam"
    bam_i = pysam.AlignmentFile(filt, "rb")
    bam_o = pysam.AlignmentFile(sin, "wb", template=bam_i)
    for r in bam_i.fetch(until_eof=True):
        if r.query_name in all_sins and r.query_name not in all_dups:
            bam_o.write(r)
    bam_i.close()
    bam_o.close()

    # 9) remap consensus
    sam_c = f"{tmp_prefix}.consensus.sam"
    rc, out, _ = sys_run(f"{options.minimap2} -ax splice {db} {fa_cons}")
    with open(sam_c, 'wb') as fw:
        fw.write(out)
    cu_uns = f"{tmp_prefix}.consensus.unsorted.bam"
    sys_run(f"{options.samtools} view -Sb {sam_c} -o {cu_uns}")
    cu_sorted = f"{tmp_prefix}.consensus.curated.bam"
    sys_run(f"{options.samtools} sort {cu_uns} -o {cu_sorted}")
    sys_run(f"{options.samtools} index {cu_sorted}")

    # 10) merge consensus + singletons → final
    final_uns = f"{tmp_prefix}.curated.unsorted.bam"
    sys_run(f"{options.samtools} cat {cu_sorted} {sin} -o {final_uns}")
    final = os.path.join(options.out_dir, f"{cell}.curated.bam")
    sys_run(f"{options.samtools} sort {final_uns} -o {final}")
    sys_run(f"{options.samtools} index {final}")

    # 11) cleanup
    if not options.keep_meta:
        for ext in (".unsorted.bam", ".filtered.bam", ".dropped.bam",
                    ".singleton.bam", ".consensus.sam", ".consensus.fa",
                    ".consensus.unsorted.bam", ".consensus.curated.bam",
                    ".curated.unsorted.bam"):
            p = tmp_prefix + ext
            if os.path.exists(p):
                os.remove(p)

    return cell

# -----------------------------------------------------------------------------
@contextmanager
def poolcontext(*args, **kwargs):
    pool = mp.Pool(*args, **kwargs)
    yield pool
    pool.terminate()

def main():
    p = argparse.ArgumentParser(
        description="Curator: map→filter→UMI-collapse→remap→final BAM"
    )
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument('--ref_genome', help='FASTA reference')
    grp.add_argument('--idx_genome', help='Minimap2 .mmi index')

    p.add_argument('--fq_dir',   required=True,
                   help='Directory of per-cell FASTQs')
    p.add_argument('--out_dir',  default='curator_res',
                   help='Where to write final BAMs')
    p.add_argument('--tmp_dir',  default='tmp', help='Scratch folder')
    p.add_argument('-t','--threads', type=int, default=1)
    p.add_argument('--umi_ld', type=int, default=1,
                   help='LD threshold for UMI merge')
    p.add_argument('--coord_buffer', type=int, default=5,
                   help='bp buffer for mapping drift')
    p.add_argument('--max_umi_duplicates', type=int, default=500)
    p.add_argument('--softclipping_thr', type=float, default=0.8)
    p.add_argument('--minimap2', default='minimap2')
    p.add_argument('--samtools', default='samtools')
    p.add_argument('--spoa',     default='spoa',
                   help='Full path to SPOA binary')
    p.add_argument('--keep_meta', action='store_true')

    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.tmp_dir, exist_ok=True)

    fqs = sorted(
        os.path.join(args.fq_dir, f)
        for f in os.listdir(args.fq_dir)
        if f.endswith('.fastq') or f.endswith('.fastq.gz')
    )
    print(f"Found {len(fqs)} FASTQs → processing with {args.threads} threads")

    start = time.time()
    with poolcontext(processes=args.threads) as pool:
        for cell in pool.map(partial(curate_one, options=args), fqs):
            if cell:
                print(f"✔️  Curated {cell}")
    dt = time.time() - start
    print(f"\nDone in {int(dt//3600)}h {int(dt%3600//60)}m {dt%60:.1f}s")

if __name__ == '__main__':
    main()
