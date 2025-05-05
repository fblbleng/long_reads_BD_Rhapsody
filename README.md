## 1. High-quality extraction of true mRNA molecules from long-read single-cell

## 🔍 Purpose
Preprocess long-read single-cell RNA-seq data to extract:

- **Cell Barcodes (CB)**
- **Unique Molecular Identifiers (UMIs)**
- **Valid insert sequences**

... from Nanopore-style reads using BD Rhapsody structure and polyA tail confirmation.

---

## 🛠️ Main Steps Performed

### 1. Adapter Detection
- Scans the **first 150 bp (head)** for adapter in **forward orientation**
- Scans the **last 150 bp (tail)** for adapter in **reverse-complemented orientation**
- Allows up to `--max_mismatch` (default: 4) using Levenshtein distance

### 2. Barcode & UMI Extraction
- Once adapter is found:
  - Extracts `--barcode_len` (default: 38 nt) immediately after adapter
  - Then extracts `--umi_len` (default: 8 nt) after that

### 3. PolyA/PolyT Tail Validation
- After CB + UMI, checks **within the next 50 nt** for a **polyA or polyT tail**
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
- Reads discarded for:
  - Adapter not found
  - Failing minimum insert length
  - Missing polyA/T
  - Other QC filters

---

## Example Read Structure Expected
```
[ADAPTER] [CB (38 nt)] [UMI (8 nt)] [INSERT] [polyA/T]
```

---

*This pipeline ensures high-quality extraction of true mRNA molecules from long-read single-cell datasets, crucial for downstream transcriptomic analyses.*

```
python3 bd_rhapsody_barcode_extractor.py -i C10_010425/subset100K.fastq.gz -o processed.fastq.gz -b barcode_list.tsv.gz -a ACACGACGCTCTTCCGATCT --barcode_len 38 --umi_len 8  --scan_region 150 --min_read_length 200 --ncores 8 
```
## 2. BD Assigner Pipeline

## 🔬 Purpose

This pipeline is designed to:
- Estimate the number of true cell barcodes (CBs) from a long-read scRNA-seq dataset
- Merge similar CBs based on sequence similarity (Levenshtein distance)
- Prepare an output table with cleaned, representative CBs

It **does NOT** collapse UMIs at this stage (UMI collapsing is to be done later).

---


## 📋 Preliminary: Generate UMI‐count Table

If you haven’t already built the CB × UMI count table from your raw `barcode_list.tsv.gz`, run:

```bash
printf "BC\tUMI\n" && \
zcat barcode_list.tsv.gz \
  | cut -f4,5 \               # extract BC and UMI
  | sort -u \                 # unique BC-UMI pairs
  | cut -f1 \                 # keep only BC
  | sort \                    # sort barcodes
  | uniq -c \                 # count UMIs per BC
  | sort -k1,1nr \            # descending by UMI count
  | awk '{print $2 "\t"$1}'   # swap and format as "BC<TAB>UMI"
| gzip > CB_count_table.tsv.gz
```
## How It Works
### 1. Load & Filter
- **Input**: the UMI-count table (`CB_count_table.tsv.gz`) with columns `BC` (barcode) and `UMI`.  
- **Filter**: drop any barcode with fewer than `--min_umi_per_cb` UMIs.  
- **Sort**: descending by UMI count.

### 2. Detect Knee & Rescue
- Compute **log₁₀(UMI + 1)** versus barcode rank.  
- Calculate per-point slope via numerical gradient.  
- **Smooth** slopes by binning on log₁₀(rank) (`--smooth_res`) and take the bin with the lowest median slope → **knee**.  
- Extend cutoff by `1 + rescue_frac` (e.g. +10 %) → **rescue**.

### 3. Select Barcodes
- Mark all barcodes from rank 1 up through the **rescue** rank as **selected_raw**.

### 4. Merge Barcodes
- Take the **selected_raw** barcodes and cluster them by Levenshtein distance:
  - **Absolute** cutoff: `--merge_dist` (e.g. ≤ 2)
  - **Fractional** cutoff: `--merge_frac` × (BC length), if specified
- Collapse each cluster into a single **representative** barcode (`rep_id`).

### 5. Output
- **Extended count table** (`--count_out`):
  - Columns:  
    - `BC`   (original barcode)  
    - `UMI`  (UMI count)  
    - `rep_id` (collapsed representative barcode)  
    - `selected_raw` (True if ≤ rescue rank)  
- **Merged list** (`--merged_out`):
  - One line per `rep_id`, listing its member barcodes.  
- **Optional plot** (`--plot`):
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



