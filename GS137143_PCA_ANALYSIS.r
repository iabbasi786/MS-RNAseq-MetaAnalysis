# ============================================================
# MS RNA-seq Meta-Analysis
# Script 1: GSE137143 - PCA Analysis (CD4+ T Cells)
# Dataset: Kim et al. 2021 (Brain) - 144 samples
#          122 MS patients + 22 Healthy Controls
#          CD4+ T cells, CD8+ T cells, CD14+ Monocytes
# Author: Iman Abbasi
# ============================================================


# ============================================================
# STEP 0: LOAD LIBRARIES
# ============================================================

library(DESeq2)
library(ggplot2)

if(!require("GEOquery")) {
  BiocManager::install("GEOquery")
  library(GEOquery)
}


# ============================================================
# STEP 1: CHECK WORKING DIRECTORY
# ============================================================
# Your working directory should be the GSE137143 data folder
# e.g. /Users/imanabbasi/Documents/University/MSc Bioinformatics/PROJECT/Data/GSE137143

cat("Working directory:", getwd(), "\n")

# Create output folders if they don't exist
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("metadata",        recursive = TRUE, showWarnings = FALSE)


# ============================================================
# STEP 2: LOAD THE COUNT MATRIX
# ============================================================
# Loads the raw gene count table
# Rows = genes (39,376 total)
# Columns = samples (419 total across all cell types)
# Values = number of RNA reads detected per gene per sample

counts_all <- read.table(
  file         = "GSE137143_raw_counts_GRCh38.p13_NCBI.tsv",
  header       = TRUE,   # First row = sample names
  row.names    = 1,      # First column = gene IDs (used as row names)
  sep          = "\t",   # Tab-separated file
  check.names  = FALSE,  # Keep sample names exactly as they are
  comment.char = ""      # Don't skip any lines
)

cat("\n=== DATA LOADED ===\n")
cat("Total genes:  ", nrow(counts_all), "\n")
cat("Total samples:", ncol(counts_all), "\n")
cat("\nFirst 5 genes, first 3 samples:\n")
print(counts_all[1:5, 1:3])


# ============================================================
# STEP 3: DOWNLOAD SAMPLE METADATA FROM GEO
# ============================================================
# Downloads sample information from the GEO database
# This tells us: which samples are MS vs HC, which cell type, etc.

cat("\n=== DOWNLOADING SAMPLE INFO FROM GEO ===\n")
cat("(This may take 1-2 minutes...)\n")

gse   <- getGEO("GSE137143", GSEMatrix = TRUE)
pheno <- pData(gse[[1]])

# Build a simple metadata data frame with Sample ID and title
metadata_all <- data.frame(
  SampleID         = rownames(pheno),
  Title            = pheno$title,
  stringsAsFactors = FALSE
)

cat("✓ Sample info downloaded\n")
cat("\nFirst 5 samples:\n")
print(head(metadata_all, 5))


# ============================================================
# STEP 4: PARSE SAMPLE INFORMATION
# ============================================================
# The GEO titles contain all the information we need, e.g.:
#   "Healthy controls, 13311d-CD4"
#   "Treatment naïve MS patients, 43213b-CD8"
# We extract Group (MS/HC) and CellType (CD4/CD8/CD14) from these

cat("\n=== PARSING SAMPLE INFORMATION ===\n")

# --- Group: MS or HC ---
# If title contains "Healthy" or "control" -> HC, otherwise -> MS
metadata_all$Group <- ifelse(
  grepl("Healthy|control", metadata_all$Title, ignore.case = TRUE),
  "HC", "MS"
)

# --- Cell Type: CD4, CD8, or CD14 ---
metadata_all$CellType <- NA

# CD4 samples: must contain CD4 but NOT CD14
# (CD14 also contains the string "CD4", so we exclude those)
metadata_all$CellType[grepl("CD4",  metadata_all$Title) &
                        !grepl("CD14", metadata_all$Title)] <- "CD4"

# CD8 samples
metadata_all$CellType[grepl("CD8",  metadata_all$Title)] <- "CD8"

# CD14 samples
metadata_all$CellType[grepl("CD14", metadata_all$Title)] <- "CD14"

# Add dataset and platform info for later use in meta-analysis
metadata_all$Dataset  <- "GSE137143"
metadata_all$Platform <- "NovaSeq6000"

cat("✓ Sample information parsed\n")
cat("\nAll samples by cell type and group:\n")
print(table(metadata_all$CellType, metadata_all$Group, useNA = "ifany"))


# ============================================================
# STEP 5: FILTER FOR CD4+ T CELLS ONLY
# ============================================================
# We focus on CD4+ T cells because:
# 1. They are central to MS immunopathology
# 2. Reduces complexity for initial analysis
# Expected: ~141 samples (120 MS + 21 HC)

cat("\n=== FILTERING FOR CD4+ T CELLS ===\n")

# Keep only CD4 rows
metadata_cd4 <- metadata_all[metadata_all$CellType == "CD4" &
                               !is.na(metadata_all$CellType), ]

# Keep only samples that exist in the count matrix
metadata_cd4 <- metadata_cd4[metadata_cd4$SampleID %in% colnames(counts_all), ]

cat("CD4+ samples found:", nrow(metadata_cd4), "\n")
cat("Group breakdown:\n")
print(table(metadata_cd4$Group))

# Subset the count matrix using match() to select only CD4 columns
# match() returns the column positions of our CD4 samples in counts_all
cd4_cols   <- match(metadata_cd4$SampleID, colnames(counts_all))
counts_cd4 <- counts_all[, cd4_cols]

# Reorder metadata rows to match the column order of the count matrix
# This is critical - rows of metadata must correspond to columns of counts
metadata_cd4 <- metadata_cd4[match(colnames(counts_cd4),
                                   metadata_cd4$SampleID), ]

cat("✓ CD4 dataset created\n")
cat("  Genes:  ", nrow(counts_cd4), "\n")
cat("  Samples:", ncol(counts_cd4), "\n")


# ============================================================
# STEP 6: QUALITY CONTROL - FILTER LOW-EXPRESSED GENES
# ============================================================
# Remove genes with very low counts - these are noise, not signal
# Rule: keep genes with >= 10 counts in at least 10 samples
# This removes genes that are barely detected across the dataset

cat("\n=== QUALITY CONTROL ===\n")
cat("Genes before filtering:", nrow(counts_cd4), "\n")

keep            <- rowSums(counts_cd4 >= 10) >= 10
counts_cd4_filt <- counts_cd4[keep, ]

cat("Genes after filtering: ", nrow(counts_cd4_filt), "\n")
cat("Genes removed:         ", nrow(counts_cd4) - nrow(counts_cd4_filt), "\n")



# ============================================================
# STEP 7: CREATE DESeq2 OBJECT
# ============================================================
# DESeq2 requires a special object that combines:
#   - The count matrix (genes x samples)
#   - The metadata (sample information)
#   - The design formula (what comparison we're making)
# design = ~ Group means we are comparing MS vs HC

cat("\n=== CREATING DESeq2 OBJECT ===\n")

# Critical safety check: metadata rows must match count matrix columns
# If this fails, something is out of order and results would be wrong
stopifnot(all(metadata_cd4$SampleID == colnames(counts_cd4_filt)))

# Make Group a factor with HC as reference level
# This means MS will be compared against HC in any downstream analysis
metadata_cd4$Group <- factor(metadata_cd4$Group, levels = c("HC", "MS"))

dds <- DESeqDataSetFromMatrix(
  countData = counts_cd4_filt,
  colData   = metadata_cd4,
  design    = ~ Group
)

cat("✓ DESeq2 object created\n")
cat("  Genes:  ", nrow(dds), "\n")
cat("  Samples:", ncol(dds), "\n")


# ============================================================
# STEP 8: NORMALISE DATA USING VST
# ============================================================
# Raw counts cannot be directly compared between samples because:
#   - Different samples have different total read counts (library size)
#   - Gene expression variance scales with mean in count data
# VST (Variance Stabilising Transformation):
#   - Corrects for library size differences
#   - Stabilises variance across expression levels
#   - Makes data suitable for PCA and visualisation
# blind = TRUE means the transformation ignores group labels
# (appropriate for QC/exploratory analysis)

cat("\n=== NORMALISING DATA (VST) ===\n")
cat("Running VST (takes 1-2 minutes)...\n")

vsd <- vst(dds, blind = TRUE)

cat("✓ Normalisation complete\n")
cat("\nExample - first gene, first 5 samples:\n")
cat("Raw counts:", as.numeric(counts(dds)[1, 1:5]), "\n")
cat("Normalised:", round(assay(vsd)[1, 1:5], 2), "\n")


# ============================================================
# STEP 9: PCA PLOT - ALL SAMPLES (INCLUDING OUTLIERS)
# ============================================================
# PCA reduces 22,000+ genes to 2 dimensions for visualisation
# PC1 = direction of greatest variation in the data
# PC2 = second greatest variation, perpendicular to PC1
# Each dot = one sample; colour = MS or HC

cat("\n=== CREATING PCA PLOT (ALL SAMPLES) ===\n")

pca_plot_all <- plotPCA(vsd, intgroup = "Group") +
  ggtitle("GSE137143: MS vs HC - CD4+ T Cells (All Samples)") +
  theme_bw() +
  theme(
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11),
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10),
    panel.grid      = element_blank()
  ) +
  scale_color_manual(
    values = c("HC" = "#3498DB", "MS" = "#E74C3C"),
    name   = "Group",
    labels = c(
      "HC" = paste0("Healthy Controls (n=", sum(metadata_cd4$Group == "HC"), ")"),
      "MS" = paste0("MS Patients (n=",      sum(metadata_cd4$Group == "MS"), ")")
    )
  )

print(pca_plot_all)

ggsave(
  filename = "results/figures/PCA_GSE137143_CD4_all.png",
  plot     = pca_plot_all,
  width    = 7, height = 6, dpi = 300, bg = "white"
)
cat("✓ Saved: results/figures/PCA_GSE137143_CD4_all.png\n")


# ============================================================
# STEP 10: IDENTIFY AND REMOVE OUTLIERS
# ============================================================
# The initial PCA shows 2 extreme outlier samples pulling PC1
# These samples have very different expression profiles to all others
# They must be removed before batch effect assessment

cat("\n=== IDENTIFYING OUTLIERS ===\n")

# Get PCA coordinates for all samples
pca_data <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)

# Samples with PC1 > 30 are clear outliers (far from the main cluster)
outliers <- pca_data[pca_data$PC1 > 30, ]
cat("Outlier samples identified:\n")
print(outliers[, c("name", "group", "PC1", "PC2")])

# Keep all samples with PC1 <= 30
keep_samples <- rownames(pca_data[pca_data$PC1 <= 30, ])
vsd_clean      <- vsd[, keep_samples]
metadata_clean <- metadata_cd4[metadata_cd4$SampleID %in% keep_samples, ]

cat("\nSamples after outlier removal:", ncol(vsd_clean), "\n")
cat("Group breakdown:\n")
print(table(metadata_clean$Group))


# ============================================================
# STEP 11: PCA PLOT - OUTLIERS REMOVED
# ============================================================
# After removing outliers, a batch effect becomes visible:
# Samples split into two clusters along PC1
# Both clusters contain MS and HC samples
# This indicates a technical (batch) effect, not biological separation

cat("\n=== CREATING PCA PLOT (OUTLIERS REMOVED) ===\n")

pca_plot_clean <- plotPCA(vsd_clean, intgroup = "Group") +
  ggtitle("GSE137143: MS vs HC - CD4+ T Cells (Outliers Removed)") +
  theme_bw() +
  theme(
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11),
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10),
    panel.grid      = element_blank()
  ) +
  scale_color_manual(
    values = c("HC" = "#3498DB", "MS" = "#E74C3C"),
    name   = "Group",
    labels = c(
      "HC" = paste0("Healthy Controls (n=", sum(metadata_clean$Group == "HC"), ")"),
      "MS" = paste0("MS Patients (n=",      sum(metadata_clean$Group == "MS"), ")")
    )
  )

print(pca_plot_clean)

ggsave(
  filename = "results/figures/PCA_GSE137143_CD4_clean.png",
  plot     = pca_plot_clean,
  width    = 7, height = 6, dpi = 300, bg = "white"
)
cat("✓ Saved: results/figures/PCA_GSE137143_CD4_clean.png\n")


# ============================================================
# STEP 12: SAVE ALL PROCESSED DATA
# ============================================================
# Save R objects so you don't need to re-run everything next time
# .rds = R data file format (like .xlsx for Excel)

cat("\n=== SAVING PROCESSED DATA ===\n")

saveRDS(counts_cd4_filt, "results/counts_cd4_filtered.rds")

saveRDS(metadata_cd4,    "metadata/metadata_GSE137143_cd4.rds")
saveRDS(dds,             "results/dds_GSE137143.rds")
saveRDS(vsd,             "results/vsd_GSE137143.rds")
saveRDS(vsd_clean,       "results/vsd_GSE137143_clean.rds")
saveRDS(metadata_clean,  "metadata/metadata_GSE137143_cd4_clean.rds")

cat("✓ All data saved\n")


# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n")
cat(rep("=", 55), "\n")

cat("              ANALYSIS COMPLETE - GSE137143\n")
cat(rep("=", 55), "\n\n")
cat("Dataset:          GSE137143 (CD4+ T cells)\n")
cat("Total samples:    ", nrow(metadata_cd4), "\n")
cat("  MS patients:    ", sum(metadata_cd4$Group == "MS"), "\n")
cat("  Healthy ctrl:   ", sum(metadata_cd4$Group == "HC"), "\n")
cat("Genes analysed:   ", nrow(dds), "\n")
cat("Outliers removed: ", nrow(pca_data) - length(keep_samples), "\n\n")
cat("Output files:\n")
cat("  results/figures/PCA_GSE137143_CD4_all.png\n")
cat("  results/figures/PCA_GSE137143_CD4_clean.png\n")
cat("  results/*.rds\n")
cat("  metadata/*.rds\n\n")
cat("Key finding: Two-cluster batch effect visible after outlier removal.\n")
cat("Both clusters contain MS and HC - separation is technical, not biological.\n")
cat("Batch correction required before integration (ComBat / limma::removeBatchEffect)\n")
