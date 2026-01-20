## TCGA-BRCA PAM50 RNA-seq Analysis
## 0. Setup
###############################################################################

setwd("C:/gdc-client_2.3_Windows_x64-py3.8-windows-2019/gdc-client_2.3_Windows_x64")
getwd()
set.seed(42)  # Set seed so results involving randomness are reproducible

# Loading required libraries
library(data.table)       # Fast reading of STAR count files
library(dplyr)            # Data manipulation and joins
library(stringr)          # Ensembl ID version stripping
library(tibble)           # Row/column name handling
library(DESeq2)           # Differential expression analysis
library(TCGAbiolinks)     # TCGA metadata (PAM50)
library(clusterProfiler)  # KEGG enrichment
library(AnnotationHub)    # Annotation resources
library(ensembldb)        # Ensembl–Entrez mapping
library(ggplot2)          # Plotting
library(ggrepel)          # Labeling volcano plots

## 1. Load cached objects if available (restart safety)

if (file.exists("rnaseq_dds_pam50.rds") &&
    file.exists("rnaseq_count_matrix_pam50_aligned.rds") &&
    file.exists("rnaseq_sample_metadata_pam50.rds") &&
    file.exists("rnaseq_gene_annotation.rds")) {
  
  dds             <- readRDS("rnaseq_dds_pam50.rds")
  count_mat       <- readRDS("rnaseq_count_matrix_pam50_aligned.rds")
  sample_meta     <- readRDS("rnaseq_sample_metadata_pam50.rds")
  gene_annotation <- readRDS("rnaseq_gene_annotation.rds")
  
} else {
  
  ###############################################################################
  ## 2. Data loading — STAR augmented gene counts
  ###############################################################################
  
  count_files <- list.files(
    path = ".",
    pattern = "rna_seq.augmented_star_gene_counts.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  stopifnot(length(count_files) > 0)
  
  sample_sheet <- fread("gdc_sample_sheet.2026-01-08.tsv")
  
  file_map <- sample_sheet %>%
    dplyr::select(file_name = `File Name`,
           sample_id = `Sample ID`,
           case_id   = `Case ID`) %>%
    inner_join(
      data.frame(file_path = count_files,
                 file_name = basename(count_files)),
      by = "file_name"
    )
  
  ###############################################################################
  ## 3. Read unstranded counts and build count matrix
  ###############################################################################
  
  read_star_unstranded <- function(file) {
    
    dt <- data.table::fread(file)
    
    dt %>%
      dplyr::filter(!grepl("^__", gene_id)) %>%
      dplyr::transmute(
        gene_id   = stringr::str_remove(gene_id, "\\.[0-9]+$"),
        gene_name = gene_name,
        gene_type = gene_type,
        count     = unstranded
      )
  }
  
  count_list <- lapply(seq_len(nrow(file_map)), function(i) {
    read_star_unstranded(file_map$file_path[i]) %>%
      dplyr::rename(!!file_map$sample_id[i] := count)
  })
  
  count_df <- Reduce(
    function(x, y) inner_join(x, y,
                              by = c("gene_id","gene_name","gene_type")),
    count_list
  )
  
  gene_annotation <- count_df %>%
    dplyr::select(gene_id, gene_name, gene_type)
  
  count_mat <- count_df %>%
    dplyr::select(-gene_name, -gene_type) %>%
    column_to_rownames("gene_id") %>%
    as.matrix()
  
  ###############################################################################
  ## 4. Metadata integration (PAM50)
  ###############################################################################
  
  pam50 <- TCGAquery_subtype("BRCA") %>%
    dplyr::select(case_id = patient,
           pam50_subtype = BRCA_Subtype_PAM50)
  
  sample_meta <- file_map %>%
    dplyr::left_join(pam50, by = "case_id") %>%
    dplyr::filter(pam50_subtype %in% c("LumA","LumB","Basal"))
  
  
  ###############################################################################
  ## 5. Sample balancing (documented step)
  ## Randomly select up to 30 samples per PAM50 subtype
  ###############################################################################
  
  sample_meta <- sample_meta %>%
    dplyr::group_by(pam50_subtype) %>%
    dplyr::group_modify(~ {
      n_take <- min(30, nrow(.x))
      dplyr::slice_sample(.x, n = n_take)
    }) %>%
    dplyr::ungroup()
  
  sample_meta$pam50_subtype <- factor(
    sample_meta$pam50_subtype,
    levels = c("LumA","LumB","Basal")
  )
  
  count_mat <- count_mat[, sample_meta$sample_id]
  dim(count_mat)
  
  ###############################################################################
  ## 6. Differential expression (DESeq2)
  ###############################################################################
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = sample_meta,
    design    = ~ pam50_subtype
  )
  
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)
  
  saveRDS(dds, "rnaseq_dds_pam50.rds")
  saveRDS(count_mat, "rnaseq_count_matrix_pam50_aligned.rds")
  saveRDS(sample_meta, "rnaseq_sample_metadata_pam50.rds")
  saveRDS(gene_annotation, "rnaseq_gene_annotation.rds")
}

###############################################################################
## 7. Differential contrasts
###############################################################################

res_BL <- results(dds, contrast = c("pam50_subtype","Basal","LumA"))
res_LB <- results(dds, contrast = c("pam50_subtype","LumB","LumA"))
res_BB <- results(dds, contrast = c("pam50_subtype","Basal","LumB"))

cat("res_BL head:\n"); print(head(res_BL))

# 11. Export DEG tables (PRIMARY DELIVERABLE)
export_deg_table <- function(res, label) {
  
  deg_df <- as.data.frame(res) %>%
    rownames_to_column("ensembl_gene_id") %>%
    dplyr::left_join(
      gene_annotation %>%
        dplyr::rename(ensembl_gene_id = gene_id),
      by = "ensembl_gene_id"
    ) %>%
    arrange(padj)
  
  cat(paste0("DEG table preview: ", label, "\n"))
  print(head(deg_df))
  
  write.csv(
    deg_df,
    paste0("rnaseq_DEG_FULL_", label, ".csv"),
    row.names = FALSE
  )
  
  write.csv(
    dplyr::filter(deg_df, !is.na(padj) & padj < 0.05),
    paste0("rnaseq_DEG_SIGNIFICANT_", label, ".csv"),
    row.names = FALSE
  )
}

export_deg_table(res_BL, "Basal_vs_LumA")
export_deg_table(res_LB, "LumB_vs_LumA")
export_deg_table(res_BB, "Basal_vs_LumB")

###############################################################################
## 9. Visualizations (QC + interpretation)
###############################################################################

## WHY:
## Visual diagnostics (QC) and interpretive plots are required to validate
## RNA-seq results and to communicate subtype-specific differences.

## WHAT BREAKS IF SKIPPED:
## Results cannot be visually verified; documentation–code mismatch.

library(ggplot2)
library(ggrepel)

###############################################################################
## 1. Variance-stabilized transformation (for visualization only)
###############################################################################

vsd <- vst(dds, blind = FALSE)

###############################################################################
## 2. PCA plot (QC across all samples)
###############################################################################

pdf("rnaseq_PCA_VST_AllSamples.pdf", width = 7, height = 6)
print(
  plotPCA(vsd, intgroup = "pam50_subtype") +
    ggtitle("PCA of TCGA-BRCA RNA-seq (VST)") +
    theme_bw()
)
dev.off()

###############################################################################
## 3. Global expression QC boxplot
###############################################################################

vsd_mat <- assay(vsd)

qc_df <- data.frame(
  expression = as.vector(vsd_mat),
  sample     = rep(colnames(vsd_mat), each = nrow(vsd_mat))
) %>%
  dplyr::left_join(
    colData(vsd) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("sample") %>%
      dplyr::select(sample, pam50_subtype),
    by = "sample"
  )

pdf("rnaseq_GlobalExpression_Boxplot_QC.pdf", width = 9, height = 5)
ggplot(qc_df, aes(x = pam50_subtype, y = expression, fill = pam50_subtype)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_bw() +
  labs(
    title = "Global Expression Distribution (VST)",
    x = "PAM50 subtype",
    y = "VST expression"
  )
dev.off()

###############################################################################
## 4. MA plots (all contrasts)
###############################################################################

pdf("rnaseq_MA_Plots_AllContrasts.pdf", width = 7, height = 6)

plotMA(res_BL, main = "MA Plot: Basal vs Luminal A")
plotMA(res_LB, main = "MA Plot: Luminal B vs Luminal A")
plotMA(res_BB, main = "MA Plot: Basal vs Luminal B")

dev.off()

###############################################################################
## 5. Volcano plot helper function (subtype-aware direction)
###############################################################################

plot_volcano <- function(res, label_left, label_right, file_name,
                         label_genes = FALSE) {
  
  df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(gene_annotation, by = "gene_id") %>%
    dplyr::mutate(
      significance = case_when(
        padj < 0.05 & log2FoldChange >= 1  ~ paste0("Up in ", label_left),
        padj < 0.05 & log2FoldChange <= -1 ~ paste0("Up in ", label_right),
        TRUE                               ~ "Not significant"
      )
    )
  
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj),
                      color = significance)) +
    geom_point(alpha = 0.6, size = 1) +
    theme_bw() +
    labs(
      title = paste("Volcano plot:", label_left, "vs", label_right),
      x = "log2 fold change",
      y = "-log10 adjusted p-value"
    ) +
    color_map <- setNames(
      c("#D55E00", "#0072B2", "grey70"),
      c(
        paste0("Up in ", label_left),
        paste0("Up in ", label_right),
        "Not significant"
      )
    )
  
  scale_color_manual(values = color_map)
  
  
  if (label_genes) {
    p <- p +
      geom_text_repel(
        data = dplyr::filter(df, padj < 0.01 & abs(log2FoldChange) >= 2),
        aes(label = gene_name),
        size = 3,
        max.overlaps = 15
      )
  }
  
  pdf(file_name, width = 7, height = 6)
  print(p)
  dev.off()
}

###############################################################################
## 5. Volcano plot helper (FIXED)
###############################################################################

plot_volcano <- function(res, label_left, label_right, file_name,
                         label_genes = FALSE) {
  
  df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(gene_annotation, by = "gene_id") %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        padj < 0.05 & log2FoldChange >= 1  ~ paste0("Up in ", label_left),
        padj < 0.05 & log2FoldChange <= -1 ~ paste0("Up in ", label_right),
        TRUE                               ~ "Not significant"
      )
    )
  
  ## Build color map safely
  color_map <- setNames(
    c("#D55E00", "#0072B2", "grey70"),
    c(
      paste0("Up in ", label_left),
      paste0("Up in ", label_right),
      "Not significant"
    )
  )
  
  scale_color_manual(values = color_map)
  
  
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj),
                      color = significance)) +
    geom_point(alpha = 0.6, size = 1) +
    theme_bw() +
    labs(
      title = paste("Volcano plot:", label_left, "vs", label_right),
      x = "log2 fold change",
      y = "-log10 adjusted p-value"
    ) +
    scale_color_manual(values = color_map)
  
  if (label_genes) {
    p <- p +
      geom_text_repel(
        data = dplyr::filter(df, padj < 0.01 & abs(log2FoldChange) >= 2),
        aes(label = gene_name),
        size = 3,
        max.overlaps = 15
      )
  }
  
  pdf(file_name, width = 7, height = 6)
  print(p)
  dev.off()
}

###############################################################################
## 6. Volcano plots (unlabeled)
###############################################################################

plot_volcano(
  res_BL, "Basal", "Luminal A",
  "rnaseq_Volcano_Unlabeled_Basal_vs_LumA.pdf"
)

plot_volcano(
  res_LB, "Luminal B", "Luminal A",
  "rnaseq_Volcano_Unlabeled_LumB_vs_LumA.pdf"
)

plot_volcano(
  res_BB, "Basal", "Luminal B",
  "rnaseq_Volcano_Unlabeled_Basal_vs_LumB.pdf"
)

###############################################################################
## 7. Labeled volcano (Basal vs Luminal A only)
###############################################################################

plot_volcano(
  res_BL, "Basal", "Luminal A",
  "rnaseq_Volcano_Labeled_Basal_vs_LumA.pdf",
  label_genes = TRUE
)

###############################################################################
## 8. DEG count bar plot
###############################################################################

###############################################################################
## DEG count bar plot (counts annotated on bars)
###############################################################################

deg_counts <- data.frame(
  Contrast = c("Basal vs LumA", "LumB vs LumA", "Basal vs LumB"),
  DEG_Count = c(
    sum(res_BL$padj < 0.05 & abs(res_BL$log2FoldChange) >= 1, na.rm = TRUE),
    sum(res_LB$padj < 0.05 & abs(res_LB$log2FoldChange) >= 1, na.rm = TRUE),
    sum(res_BB$padj < 0.05 & abs(res_BB$log2FoldChange) >= 1, na.rm = TRUE)
  )
)

pdf("rnaseq_DEG_Counts_Barplot.pdf", width = 6, height = 4)

ggplot(deg_counts, aes(x = Contrast, y = DEG_Count, fill = Contrast)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = DEG_Count),
    vjust = -0.4,
    size = 4
  ) +
  theme_bw() +
  labs(
    title = "Number of Significant DEGs per Contrast",
    y = "DEG count (padj < 0.05, |log2FC| ≥ 1)",
    x = NULL
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  expand_limits(y = max(deg_counts$DEG_Count) * 1.1)

dev.off()

###############################################################################
## KEGG pathway enrichment (Ensembl → Entrez, namespace-safe)
###############################################################################

## WHY:
## KEGG requires Entrez Gene IDs. TCGA STAR outputs are Ensembl-based.
## Explicit mapping avoids silent ID errors and reviewer concerns.

## WHAT BREAKS IF SKIPPED:
## KEGG will fail or return empty results due to unmapped identifiers.

library(AnnotationHub)
library(ensembldb)
library(clusterProfiler)

## ---------------------------------------------------------------------------
## 1. Retrieve Ensembl–Entrez mapping from AnnotationHub
## ---------------------------------------------------------------------------

ah  <- AnnotationHub()
edb <- ah[["AH73881"]]  # Ensembl 97 EnsDb for Homo sapiens

ens_entrez <- genes(
  edb,
  columns = c("gene_id", "entrezid"),
  filter  = GeneIdFilter(keys(edb))
) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, entrezid) %>%
  dplyr::filter(!is.na(entrezid)) %>%
  dplyr::distinct()

## Sanity checks
cat("Ensembl–Entrez mapping dim:\n")
print(dim(ens_entrez))
print(head(ens_entrez))

## ---------------------------------------------------------------------------
## 2. Define gene universe (explicit, consistent across contrasts)
## ---------------------------------------------------------------------------

gene_universe <- rownames(dds)

universe_entrez <- data.frame(gene_id = gene_universe) %>%
  dplyr::inner_join(ens_entrez, by = "gene_id") %>%
  dplyr::pull(entrezid) %>%
  unique()

cat("Gene universe (Entrez) size:\n")
print(length(universe_entrez))

## ---------------------------------------------------------------------------
## 3. KEGG enrichment function (namespace-safe)
## ---------------------------------------------------------------------------

run_kegg <- function(res, label) {
  
  ## Extract significant DEGs (documented threshold)
  deg_sig <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::filter(
      !is.na(padj),
      padj < 0.05
    )
  
  ## Map DEGs to Entrez IDs
  deg_entrez <- deg_sig %>%
    dplyr::inner_join(ens_entrez, by = "gene_id") %>%
    dplyr::pull(entrezid) %>%
    unique()
  
  ## Defensive check
  if (length(deg_entrez) == 0) {
    message("No Entrez-mapped DEGs for ", label, ". KEGG skipped.")
    return(NULL)
  }
  
  ## Run KEGG enrichment
  ek <- enrichKEGG(
    gene         = deg_entrez,
    universe     = universe_entrez,
    organism     = "hsa",
    pvalueCutoff = 0.05
  )
  
  ## Save results
  write.csv(
    as.data.frame(ek),
    paste0("rnaseq_KEGG_", label, ".csv"),
    row.names = FALSE
  )
  
  return(ek)
}

## ---------------------------------------------------------------------------
## 4. Run KEGG for all documented contrasts
## ---------------------------------------------------------------------------

kegg_BL <- run_kegg(res_BL, "Basal_vs_LumA")
kegg_LB <- run_kegg(res_LB, "LumB_vs_LumA")
kegg_BB <- run_kegg(res_BB, "Basal_vs_LumB")



###############################################################################
## END OF PIPELINE
###############################################################################
