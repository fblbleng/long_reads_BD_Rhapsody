### 1. High-quality extraction of true mRNA molecules from long-read single-cell

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
# 2. BD Rhapsody-Style Long-Read CB Assigner

## 🔍 Purpose
This pipeline refines the list of extracted cell barcodes (CBs) and unique molecular identifiers (UMIs) from long-read single-cell RNA sequencing (scRNA-seq) datasets (e.g., Nanopore) using BD Rhapsody protocols.

It:
- Estimates the number of true cells.
- Collapses similar CBs with small sequencing errors.
- Outputs cleaned barcode tables for downstream analysis.

---

## 🛠️ Major Steps

### 1. Load Barcode/UMI List
- Input a TSV or TSV.GZ file containing fields like `read_id`, `orientation`, `BC`, `UMI`, etc.

### 2. UMI Counting
- Counts the number of **unique UMIs** per CB.

### 3. Cell Number Estimation
- Estimates the number of real cells via "knee" detection on the log-scaled UMI count distribution.
- Optionally produces a plot of UMI counts vs barcode rank (`--plot`).

### 4. Barcode Merging
- Collapses barcodes within a Levenshtein distance (default: 3 mismatches) to rescue noisy or erroneous CBs (`--merge_distance`).

### 5. Filter Cells
- Retains CBs with at least `--min_UMI_per_CB` UMIs (default: 2).
- (Planned) Filters based on number of genes detected per CB (`--min_genes`).

### 6. Output
- Produces a clean `.tsv.gz` table associating each read to its merged "representative" CB.

---

## 📂 Output Files
- **CB_count_table.tsv.gz**: cleaned barcode assignment with columns:
  - read_id
  - orientation
  - BC_start
  - BC (original)
  - UMI
  - Seq_end
  - mean_BC_quality
  - Representative_CB (merged barcode)

- **knee_plot.png** (optional): Visualizes the barcode "knee" for cell number estimation.

---

## 🚀 Example Command
```bash
python3 assigner.py \
  -i barcode_list.tsv.gz \
  -o CB_count_table.tsv.gz \
  --min_UMI_per_CB 2 \
  --merge_distance 3 \
  --plot \
  -t 4
```

---

## 📊 Parameters Summary
| Parameter              | Description                                           | Default |
|------------------------|-------------------------------------------------------|---------|
| `-i`                   | Input barcode file                                    | Required |
| `-o`                   | Output file                                           | CB_count_table.tsv.gz |
| `--min_UMI_per_CB`      | Min UMIs to keep a CB                                 | 2 |
| `--min_genes`           | (Reserved) Min genes to keep a CB (not yet applied)   | 300 |
| `--merge_distance`      | Levenshtein distance for merging barcodes             | 3 |
| `--dynamic_rescue`      | (Reserved) Rescue true cells beyond knee             | Off |
| `--plot`                | Output knee plot                                      | Off |
| `-t, --threads`         | Threads for parallelization                          | 4 |

---

## 📈 Notes
- **Flexible and Modular:**
  - Supports gzipped inputs and outputs.
  - Adjustable sensitivity for CB collapsing.
- **No Whitelist Needed:**
  - Completely unsupervised; does not rely on predefined barcode lists.

---

## ✨ Future Improvements
- Gene filtering based on expression matrices.
- More dynamic knee detection (with second derivative smoothing).
- UMI collapsing with distance-based correction.
- Output per-cell FASTQ deconvolution.

---

# ⚡ Thank you for using this assigner pipeline!
Feel free to suggest improvements or request new features!

