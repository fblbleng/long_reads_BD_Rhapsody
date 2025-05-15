## 1. High-quality extraction of true mRNA molecules from long-read single-cell

## 🔍 Purpose
Preprocess long-read single-cell RNA-seq data to extract:

- ###Cell Barcodes (CB)###
- ###Unique Molecular Identifiers (UMIs)###
- ###Valid insert sequences###

... from Nanopore-style reads using BD Rhapsody structure and polyA tail confirmation.

---

## 🛠️ Main Steps Performed

### 1. Adapter Detection
- Scans the **first 150 bp** (head)for adapter in **forward orientation**
- Scans the **last 150 bp** (tail) for adapter in **reverse-complemented orientation**
- Allows up to `--max_mismatch` (default: 4) using Levenshtein distance

### 2. Barcode & UMI Extraction
- Once adapter is found:
  - Extracts `--barcode_len` (default: 38 nt) immediately after adapter
  - Then extracts `--umi_len` (default: 8 nt) after that

### 3. PolyA/PolyT Tail Validation
- After CB + UMI, checks ###within the next 50 nt### for a ###polyA or polyT tail###
- Ensures the read likely came from an mRNA transcript

### 4. Insert Trimming & Filtering
- Keeps only insert sequence after adapter + CB + UMI
- Discards if insert is shorter than `--min_read_length` (default: 200 nt)
- Discards if no valid polyA/polyT tail detected after CB/UMI

### 5. Output Files
- ✅ `processed.fastq.gz`: Cleaned reads with adapter/CB/UMI removed
- 📄 `barcode_list.tsv.gz`: Per-read table with:
  - `read_id`, `orientation`, `BC_start`, `BC`, `UMI`, `Seq_end`, `mean_BC_quality`

### 6. Statistics
Printed to console after processing:
- Total reads
- Reads with valid adapter
---

## Example Read Structure Expected
```
[ADAPTER] [CB (38 nt)] [UMI (8 nt)] [INSERT] [polyA/T]
```

---

*This pipeline ensures high-quality extraction of true mRNA molecules from long-read single-cell datasets, crucial for downstream transcriptomic analyses.*

```
python3 bd_rhapsody_barcode_extractor.py -i ~/Disk/Donnees_brutes/Rhapsody_PC/LongRead/C8_C9/C9.fastq.gz  -o processed_ALL.fastq.gz -b barcode_ALL_list.tsv.gz -a ACACGACGCTCTTCCGATCT --barcode_len 35 --umi_len 8  --scan_region 150 --min_read_length 200 --ncores 40 --chunksize 10000
```

## 2. BD Assigner Pipeline
`

## 🔬 Purpose

This pipeline is designed to:
- Estimate the number of true cell barcodes (CBs) from a long-read scRNA-seq dataset
- Merge similar CBs based on sequence similarity (Levenshtein distance)
- Prepare an output table with cleaned, representative CBs

It ###does NOT### collapse UMIs at this stage (UMI collapsing is to be done later).

---

## How It Works
### 1. Load & Filter
- ###Input###: the UMI-count table (`CB_count_table.tsv.gz`) with columns `BC` (barcode) and `UMI`.  
- ###Filter###: drop any barcode with fewer than `--min_umi_per_cb` UMIs.  
- ###Sort###: descending by UMI count.

### 2. Detect Knee & Rescue
- Compute ###log₁₀(UMI + 1)### versus barcode rank.  
- Calculate per-point slope via numerical gradient.  
- ###Smooth### slopes by binning on log₁₀(rank) (`--smooth_res`) and take the bin with the lowest median slope → ###knee###.  
- Extend cutoff by `1 + rescue_frac` (e.g. +10 %) → ###rescue###.

### 3. Select Barcodes
- Mark all barcodes from rank 1 up through the ###rescue### rank as ###selected_raw###.

### 4. Merge Barcodes
- Take the ###selected_raw### barcodes and cluster them by Levenshtein distance:
  - ###Absolute### cutoff: `--merge_dist` (e.g. ≤ 2)
  - ###Fractional### cutoff: `--merge_frac` × (BC length), if specified
- Collapse each cluster into a single ###representative### barcode (`rep_id`).

### 5. Output
- ###Extended count table### (`--count_out`):
  - Columns:  
    - `BC`   (original barcode)  
    - `UMI`  (UMI count)  
    - `rep_id` (collapsed representative barcode)  
    - `selected_raw` (True if ≤ rescue rank)  
- ###Merged list### (`--merged_out`):
  - One line per `rep_id`, listing its member barcodes.  
- ###Optional plot### (`--plot`):
  - Saves `knee_rescue.png` showing the log-log UMI curve with knee & rescue cutoffs.

---

## 🔠 Usage

```bash
python3 bd_assign.py \
  -i CB_count_table.tsv.gz \
  -o CB_counting_ext.tsv.gz \
  -m CB_merged_list.txt \
  [--min_umi_per_cb 5] \
  [--smooth_res 0.001] \
  [--rescue_frac 0.10] \
  [--merge_dist 2 | --merge_frac 0.2] \
  [--forced_no 0] \
  [--plot]
```

### Arguments:
- `-i`  : Input barcode list
- `-o`  : Output CB count table (default `CB_count_table.tsv.gz`)
- `--min_UMI_per_CB` : Minimum UMIs per CB to consider (default: 2)
- `--merge_distance` : Maximum distance to merge CBs (default: 3)
- `--plot` : Whether to plot knee curve
- `-t`  : Number of threads (default: 4)

---

## 📊 Output Files

- `CB_count_table.tsv.gz`
  - Cleaned list of barcodes with a column `Representative_CB`

- `knee_plot.png`
  - Visualizes the UMI distribution and detected knee point (optional)

---

## 🚀 Notes

- Only barcodes are collapsed, no UMI collapsing yet.
- Filtering on `--min_genes` is NOT performed here; only UMI count is used.
- For better results, especially with noisy data like BD Rhapsody, you can adjust `--merge_distance` and `--min_UMI_per_CB` accordingly.

---

# 3. BD Rhapsody-Style Parallel Demultiplexer

## 🔍 Purpose
Demultiplex a cleaned, barcode‐trimmed FASTQ into per‐cell FASTQ files by matching each read’s extracted cell barcode (CB) and UMI against a set of ### representative barcodes ###.  
Uses fuzzy matching (Levenshtein distance) and parallel processing for speed.

---

## 🛠️ Main Steps

### 1. Load Representative Barcodes  
Read the list of collapsed “true” barcodes from `--rep_file` (e.g. `CB_count_table.tsv.gz`, column `Representative_CB`).

### 2. Load Read → Barcode Mapping  
Read the per‐read table `--mapping` (e.g. `barcode_list.tsv.gz`) 
Build a dictionary `read_id → (raw_CB, UMI)`.

### 3. Pre-map Raw CB → Representative CB  
In parallel (using `--threads` workers), for each `read_id`:
- Compute Levenshtein distance between its raw CB and each rep CB
- Assign the read to the rep with the smallest distance ≤ `--max_mismatch`  

This produces a `premap` lookup: only reads within tolerance are kept.

### 4. Demultiplex FASTQ  
Stream the trimmed FASTQ (`--input_fastq`) in ###chunks###:
- For each record, look up its `read.id` in `premap`
- Annotate the header:
- Write the record into `output_dir/<rep_CB>.fastq.gz`

All I/O is gzip-compatible and thread-safe.

---

## 📂 Input & Output

### Inputs
- `-i, --input_fastq`  
Cleaned FASTQ from the extraction step (e.g. `processed.fastq.gz`)
- `-m, --mapping`  
TSV(.gz) table mapping `read_id → BC, UMI` (from `barcode_list.tsv.gz`)
- `-r, --rep_file`  
TSV(.gz) table of representative barcodes (e.g. `CB_count_table.tsv.gz`)

### Outputs
- `-o, --output_dir`  
Directory of per-cell FASTQ files, one `*.fastq.gz` per rep CB
### Console summary:  
- Total reads scanned  
- Reads successfully demultiplexed  
- Number of output files (cells)

---

## ⚙️ Parameters

| Option                    | Description                                                                                     | Default |
|---------------------------|-------------------------------------------------------------------------------------------------|---------|
| `-i, --input_fastq`       | Processed FASTQ to demultiplex (gzip OK)                                                        | (none)  |
| `-m, --mapping`           | Read‐to‐barcode table (`read_id`, `BC`, `UMI`)                                                  | (none)  |
| `-r, --rep_file`          | Representative CB list (`Representative_CB` column)                                            | (none)  |
| `-o, --output_dir`        | Output folder for per-cell FASTQ files                                                          | `demux_fastq` |
| `--max_mismatch`          | Maximum Levenshtein distance for CB→rep matching                                                | `1`     |
| `-t, --threads`           | Number of parallel workers for pre-mapping and demux                                              | `4`     |

---

## 🔧 Example Usage

```bash
python3 demultiplex.py \
-i processed.fastq.gz \
-m barcode_list.tsv.gz \
-r CB_count_table.tsv.gz \
-o demuxed_fastq \
--max_mismatch 2 \
-t 8
```
### 4 Curation step  takes per-cell demultiplexed FASTQ files (one per cell) and produces ###“curated”### BAMs by:

1. *** Splice-aware alignment *** (Minimap2)  
2. *** Soft-clip filtering *** (remove poorly aligned reads)  
3. ###UMI clustering & consensus### (SPOA)  
4. ###Re-mapping consensus### (Minimap2)  
5. ###Merging consensus + singleton reads### → final sorted/indexed BAM  

---

### Reference indexing
Build a .mmi for splice mode:

```bash
minimap2 -d GRCh38.mmi Homo_sapiens.GRCh38.dna.primary_assembly.fa
```
Use the same genome (and matching GTF) downstream for quantification.

## Inputs

- ###*`--fq_dir`###  
  Directory containing per-cell FASTQ files (e.g. `demux_fastq/CTTCAG...TTACTT.fastq.gz`).

- ###*Reference### (one of):
  - ######`--ref_genome`###### — Path to genome FASTA  
  - ######`--idx_genome`###### — Path to prebuilt Minimap2 index (`.mmi`)

- ###*Directories###:
  - ###`--tmp_dir`### — Temporary working folder (default: `tmp`)  
  - ###`--out_dir`### — Final BAM output folder (default: `curator_res`)  

---

## Dependencies

- ###Python 3###  
- ###pysam### (for BAM I/O)  
- ###Biopython### (`Bio.SeqIO`)  
- ###Minimap2### (CLI in `$PATH`)  
- ###samtools### (CLI in `$PATH`)  
- ###spoa### (CLI in `$PATH`)  

---

## Installation

```bash
# Using conda for all dependencies:
conda install -c bioconda pysam biopython minimap2 samtools spoa
🚀 Usage
bash
Copy
Edit
python3 curator.py \
  --fq_dir     demux_fastq          \
  --idx_genome /path/to/GRCh38.mmi  \
  --tmp_dir    bdrhapsody_gps_res/tmp \
  --out_dir    bdrhapsody_gps_res/final_bams \
  --threads    8                     \
  --umi_ld     2                     \
  --coord_buffer 5                   \
  --max_umi_duplicates 500           \
  --softclipping_thr 0.8             \
  [--keep_meta]                      \
  [--ref_genome /path/to/GRCh38.fa]  # if no --idx_genome
--threads: number of parallel workers

--umi_ld: merge UMIs within this Levenshtein distance

--coord_buffer: bp drift buffer when clustering by mapping position

--max_umi_duplicates: cap per-UMI cluster size used to build consensus

--softclipping_thr: require (matches)/(matches+softclips) ≥ threshold

--keep_meta: retain all intermediate BAM/SAM files
```
## Step-by-Step Details

### 1. Discover FASTQs
- Scans `--fq_dir` for `*.fastq` / `*.fastq.gz` files.

### 2. Parallel per-cell processing
- Spawns `curate_one()` across cells using `multiprocessing.Pool`.

---

### 3. `curate_one(fq_path)`
  
### a. Map → Unsorted BAM ###  
```bash
minimap2 -ax splice <ref> <fastq> \
  | samtools view -Sb - > tmp/<cell>.unsorted.bam
```
### b. Soft-clip filtering ###  
- Keeps reads with a high fraction of aligned bases  
- Writes dropped reads to `tmp/<cell>.dropped.bam`  

### c. Index filtered BAM ###  
```bash
samtools index tmp/<cell>.filtered.bam
```
### d. Load raw sequences

- Builds a read_id → sequence dictionary from the FASTQ

### e. Cluster by mapping position

Groups reads on the same chromosome within --coord_buffer bp windows

###f . Collapse UMIs (per cluster)

-Extract UMI tag from read_id

-Merge UMI groups within Levenshtein distance ≤ --umi_ld

-Write a tiny FASTA and run spoa → consensus sequence

### Collect:

-consensus: polished sequences

-duplicates: reads used to build consensus

-singletons: UMIs with only one read

### g. Write consensus FASTA

tmp/<cell>.consensus.fa

###  h. Extract singleton BAM

- From tmp/<cell>.filtered.bam

### i. Remap consensus → sorted/indexed BAM

```bash

minimap2 -ax splice <ref> tmp/<cell>.consensus.fa \
  | samtools view -Sb - \
  | samtools sort -o tmp/<cell>.consensus.curated.bam
samtools index tmp/<cell>.consensus.curated.bam
```
### j. Merge consensus + singleton BAM → final sorted/indexed BAM

```bash

samtools cat tmp/<cell>.consensus.curated.bam tmp/<cell>.singleton.bam \
  | samtools sort -o <out_dir>/<cell>.curated.bam
samtools index <out_dir>/<cell>.curated.bam
```
### k. Clean up

- Removes all intermediate files unless --keep_meta is set

# Reporter Expression Pipeline

## 🔍 Purpose
This script generates a gene-by-cell expression matrix from per-cell curated BAM files produced by the Curator step. It:

- Counts reads overlapping annotated genes using **featureCounts**
- Aggregates per-cell counts into a single matrix
- Filters out low-quality cells with too few detected genes

## 📂 Input Requirements

### 1. Per-cell curated BAMs
- Directory of `*.curated.bam` files (each must be indexed with `.bai`)
- Produced by `curator.py` into `<out_dir>` (default: `bdrhapsody_gps_res`)

### 2. GTF annotation (`--gtf`)
- A standard gene annotation file (e.g., Gencode or Ensembl GTF) matching your reference genome

### 3. Filtered barcode list (`--sel_bc_o`)
- Plain-text file listing one accepted cell barcode per line
- Only BAMs whose cell ID appears in this list will be quantified

## 🛠️ Workflow

1. **Discover curated BAMs**  
   Scans `<O_DIR>` for `*.curated.bam` files.

2. **Optional barcode filtering**  
   If `--sel_bc_o` is provided, restrict to those cell IDs; otherwise process all BAMs.

3. **Run featureCounts**  
   For each cell:

   ``` bash
   featureCounts \
     -T <threads> \
     -a <gtf> \
     -o <tmp_dir>/<cell>.counts.txt \
     <out_dir>/<cell>.curated.bam

## Detailed Workflow Steps

### FeatureCounts Execution
- Uses `featureCounts` to quantify reads per gene per cell
- **Custom binary path**: Use `--featurecounts` to override the path to the featureCounts binary if not in system PATH

### Aggregate Counts
1. Parse each `<cell>.counts.txt` file
2. Build a sparse matrix with:
   - Rows: `gene_id` (from GTF)
   - Columns: Cell barcodes
   - Values: Raw read counts

### Quality Filtering
- Removes any cell (column) with fewer than `--min_gene_no` detected genes
- Filter threshold is configurable (default: 300 genes)

### Output Generation
- Writes final matrix as TSV:
  - File: `<O_DIR>/<o_name>` (default: `matrix.tsv`)
  - Format:
    ```
    gene_id    cell1    cell2    ...
    ENSG0001   15       0        ...
    ENSG0002   3        42       ...
    ```

### Logging & Metadata
- Progress and timing logged to `--log` (default: `bdrhapsody_gps_res/logs/reporter_expression.log.txt`)
- Intermediate files:
  - Kept in `--tmp_dir` if `--keep_meta` is set
  - Otherwise automatically cleaned up

## Command-Line Options

| Flag | Description | Default |
|------|-------------|---------|
| `-h`, `--help` | Show help message | |
| `-d O_DIR` | Output directory | `bdrhapsody_gps_res` |
| `--tmp_dir TMP_DIR` | Temp directory for counts | `tmp` |
| `--gtf GTF` | (Required) Gene annotation GTF | |
| `-o O_NAME` | Output matrix filename | `matrix.tsv` |
| `--log LOG_F_NAME` | Log file path | `bdrhapsody_gps_res/logs/reporter_expression.log.txt` |
| `-t NCORES` | CPU threads | `1` |
| `--min_gene_no MIN_GENE_NO` | Min genes per cell | `300` |
| `--sel_bc_o SEL_BC_O` | Barcode whitelist file | `filtered_barcode_list.txt` |
| `--keep_meta` | Keep temp count files | `False` |
| `--featurecounts FEATURECOUNTS` | Path to featureCounts binary | `featureCounts` |

## Example Usage

```bash
python3 reporter_expression.py \
  -d bdrhapsody_gps_res \
  --tmp_dir bdrhapsody_gps_res/tmp \
  --gtf /path/to/gencode.v44.annotation.gtf \
  -o gene_matrix.tsv \
  --log bdrhapsody_gps_res/logs/expression.log.txt \
  -t 8 \
  --min_gene_no 50 \
  --sel_bc_o filtered_barcodes.txt \
  --featurecounts /usr/local/bin/featureCounts
```
