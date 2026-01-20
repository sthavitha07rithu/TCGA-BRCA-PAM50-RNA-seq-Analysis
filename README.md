# TCGA-BRCA PAM50 RNA-seq Differential Expression Analysis

## Overview

This repository contains a **restart-safe, end-to-end RNA-seq analysis pipeline** for **TCGA Breast Invasive Carcinoma (TCGA-BRCA)**, focused on identifying **differentially expressed genes (DEGs) across PAM50 molecular subtypes**.

The analysis is implemented in a single R script (`mini_project_FINAL.R`) and produces **gene-level DEG tables as the primary deliverable**, with quality-control visualizations and KEGG pathway enrichment as supporting analyses.

---

## Dataset

- **Project:** TCGA-BRCA
- **Data type:** RNA-seq gene expression quantification
- **Workflow:** STAR – augmented gene counts
- **Counts used:** Unstranded gene-level counts
- **Subtypes analyzed:** PAM50  
  - Luminal A (LumA)  
  - Luminal B (LumB)  
  - Basal  

Raw data are downloaded using the **GDC Data Transfer Tool** and are not included in this repository.

---

## Analysis Workflow

The pipeline follows the workflow implemented in `mini_project_FINAL.R`:

1. **Data loading**
   - STAR augmented gene count files are discovered recursively
   - TCGA sample metadata are read from the GDC sample sheet

2. **Count matrix construction**
   - Unstranded counts are extracted
   - Ensembl gene version suffixes are removed
   - Gene annotation (Ensembl ID, gene symbol, gene type) is preserved
   - Only genes consistently present across samples are retained (conservative approach)

3. **Metadata integration**
   - PAM50 subtype annotations are retrieved using `TCGAbiolinks`
   - Samples are restricted to LumA, LumB, and Basal subtypes

4. **Sample balancing**
   - Up to **30 samples per PAM50 subtype** are randomly selected
   - The random selection is reproducible via a fixed seed

5. **Differential expression analysis**
   - Performed using **DESeq2**
   - Design formula: `~ pam50_subtype`
   - Low-count genes are filtered prior to model fitting

6. **Contrasts tested**
   - Basal vs Luminal A  
   - Luminal B vs Luminal A  
   - Basal vs Luminal B  

7. **Primary deliverables**
   - Full DEG tables for each contrast
   - Statistically significant DEG tables (adjusted p-value < 0.05)

8. **Visualizations**
   - PCA (VST-transformed counts)
   - Global expression QC boxplot
   - MA plots for all contrasts
   - Volcano plots (labeled and unlabeled)
   - DEG count bar plot with counts annotated

9. **KEGG pathway enrichment**
   - Ensembl → Entrez ID mapping via `AnnotationHub` / `ensembldb`
   - Enrichment performed using `clusterProfiler`
   - Same gene universe used across all contrasts

10. **Restart safety**
    - Intermediate objects are saved as `.rds` files
    - If present, cached objects are reloaded to avoid recomputation

---

## Primary Outputs

### Differential Expression Tables (Primary Deliverable)

For each contrast:

- `rnaseq_DEG_FULL_<contrast>.csv`  
  Full DESeq2 results with gene annotation

- `rnaseq_DEG_SIGNIFICANT_<contrast>.csv`  
  Subset of DEGs with adjusted p-value < 0.05

These files represent the **authoritative gene-level results** of the analysis.

---

### Figures

All figures are saved as PDF files prefixed with `rnaseq_`, including:

- PCA of VST-transformed counts
- Global expression QC boxplot
- MA plots
- Volcano plots (labeled and unlabeled)
- DEG count bar plot

---

### KEGG Enrichment Results

For each contrast:

- `rnaseq_KEGG_<contrast>.csv`

These tables contain pathway-level enrichment statistics and are intended as **supporting, interpretive results**, not the primary outcome.

---

## Running the Pipeline

1. Download TCGA-BRCA STAR count data using the GDC Data Transfer Tool.
2. Place the downloaded files and the GDC sample sheet in the working directory.
3. Open R and set the working directory to the project folder.
4. Run:

```r
source("mini_project_FINAL.R")
```

The script will either:
Load previously saved .rds objects, or
Rebuild the analysis from raw inputs if cached files are absent.

## Reproducibility
The pipeline is deterministic given identical inputs
Random sampling is controlled with a fixed seed
Intermediate results are saved and reloadable
Sanity checks (dim(), head()) are printed during execution
Internet access is required only for PAM50 metadata retrieval and KEGG annotation

## Software Requirements
R (≥ 4.2 recommended)
CRAN packages:
data.table
dplyr
stringr
tibble
ggplot2
ggrepel

Bioconductor packages:
DESeq2
TCGAbiolinks
clusterProfiler
AnnotationHub
ensembldb

## Data Download Instructions (GDC Data Transfer Tool)

### Prerequisites

To reproduce the analysis, RNA-seq count files must be downloaded from the **Genomic Data Commons (GDC)**.

You will need:
- A GDC account: https://portal.gdc.cancer.gov/
- The **GDC Data Transfer Tool (CLI)**
- A GDC **manifest file** (`gdc_manifest.txt`)

---

### Step 1: Install the GDC Data Transfer Tool

Download the appropriate version for your operating system from:

https://gdc.cancer.gov/access-data/gdc-data-transfer-tool

Verify installation:

```bash
gdc-client --version
```

### Step 2: Obtain the Manifest File

1. Visit the GDC Data Portal
2. Apply the following filters:
   - Project: TCGA-BRCA
   - Data Type: Gene Expression Quantification
   - Workflow Type: STAR - Counts
3. Add files to the cart
4. Download the manifest file (gdc_manifest.txt)

Place the manifest file in the root directory of this repository.

Important: Do not commit access tokens or controlled-access files to GitHub.

---

### Step 3: Download the Data via CLI

From the repository root directory, run:

gdc-client download -m gdc_manifest.txt

This will:
- Create UUID-named subdirectories
- Download augmented STAR gene count files
- Preserve the directory structure expected by the analysis script

The pipeline automatically discovers these files recursively.

---

### Step 4: Verify Download Structure

After download, the directory structure should resemble:

.
├── gdc_manifest.txt
├── <UUID_1>/
│   └── *.rna_seq.augmented_star_gene_counts.tsv
├── <UUID_2>/
│   └── *.rna_seq.augmented_star_gene_counts.tsv
└── ...

No manual renaming or file movement is required.
<img width="814" height="1568" alt="image" src="https://github.com/user-attachments/assets/4cc3ef00-d996-4973-a5ff-d908ba73489d" />

