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

## 🔮 Example Read Structure Expected
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

## 🔧 How It Works

### 1. Load CB List
- Reads the `barcode_list.tsv.gz` generated from preprocessing (Scanner pipeline)

### 2. Estimate Cell Number
- Computes UMI counts per CB
- Plots UMI count vs barcode rank (log-log)
- Detects "knee point" (where signal drops dramatically) to estimate the number of cells

### 3. Merge Barcodes
- Among top-ranked CBs (up to knee point), barcodes are merged if:
  - Levenshtein distance ≤ `--merge_distance` (default: 3)

### 4. Output
- A cleaned table with one **Representative_CB** per cell
- A plot (`knee_plot.png`) showing the UMI distribution and knee point if `--plot` is enabled

---

## 📂 Input Requirements

- `barcode_list.tsv.gz`
  - Must contain at least columns: `read_id`, `BC`, `UMI`

---

## 🔠 Usage

```bash
python3 bd_assigner.py \
    -i barcode_list.tsv.gz \
    -o CB_count_table.tsv.gz \
    --min_UMI_per_CB 10 \
    --merge_distance 3 \
    --plot \
    -t 4
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



