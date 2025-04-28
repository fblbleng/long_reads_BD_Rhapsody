# High-quality extraction of true mRNA molecules from long-read single-cell

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
