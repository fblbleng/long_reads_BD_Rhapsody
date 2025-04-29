import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === Load your CB UMI count table ===
umi_table = pd.read_csv("CB_count_table.tsv.gz", sep="\t")

# === Force UMI to numeric ===
umi_table["UMI"] = pd.to_numeric(umi_table["UMI"], errors="coerce").fillna(0).astype(int)

# === Sort barcodes by UMI counts ===
umi_table = umi_table.sort_values("UMI", ascending=False).reset_index(drop=True)

# === Compute log10(UMI) for nicer plot ===
umi_table["log10_UMI"] = np.log10(umi_table["UMI"] + 1)

# === Compute derivatives (slopes) ===
umi_table["slope"] = umi_table["log10_UMI"].diff()

# === Crude Knee Detection (sharpest drop) ===
knee_idx = umi_table["slope"].idxmin()

# === Dynamic Rescue: extend 10% more barcodes ===
rescue_idx = int(knee_idx * 1.10)

# === Plot ===
plt.figure(figsize=(10,6))
plt.plot(umi_table.index, umi_table["log10_UMI"], label="log10(UMI) per CB")
plt.axvline(knee_idx, color="red", linestyle="--", label=f"Knee Point: {knee_idx}")
plt.axvline(rescue_idx, color="green", linestyle="--", label=f"Dynamic Rescue: {rescue_idx}")
plt.xlabel("Barcode Rank")
plt.ylabel("log10(UMI Count)")
plt.title("Crude Knee vs Dynamic Rescue on CBs")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
