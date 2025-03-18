# 🧬 RNAseq Pipeline for C. elegans Stress Response

## 📌 Overview
This project analyzes RNAseq data from _C. elegans_ to identify differential gene expression but not the data.  
It includes transcript quantification, differential expression analysis, and visualization.


---

## 📂 Repository Structure
- **`Code/`** → Contains all scripts for data analysis.
  - `salmon_quant.sh` → Runs Salmon for transcript quantification.
  - `dseq_analysis.R` → Performs DESeq2 differential expression analysis.
  - `ex_heatmap_volcano.R` → Creates Volcano Plot & Heatmap.
- **`Visualizations/`** → Contains the output figures.
  - `N2_on_strep.png` → Volcano plot.
  - `stress_genes_N2_sgord.png` → Heatmap.

---
