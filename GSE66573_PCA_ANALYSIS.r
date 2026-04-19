# ============================================================
# MS RNA-seq Meta-Analysis
# Script 3: GSE66573 - PCA Analysis
# Dataset: 14 samples - MS patients vs Healthy Controls
#          Whole blood RNA-seq
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
# Your working directory should be the GSE66573 data folder

cat("Working directory:", getwd(), "\n")

# Create output folders
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("metadata",        recursive = TRUE, showWarnings = FALSE)


# ============================================================
# STEP 2: LOAD THE COUNT MATRIX
# ============================================================
# Small dataset: 14 samples

counts_all <- read.table(
  file         = "GSE66573_raw_counts_GRCh38.p13_NCBI.tsv",
  header       = TRUE,
  row.names    = 1,
  sep          = "\t",
  check.names  = FALSE,
  comment.char = ""
)

cat("\n=== DATA LOADED ===\n")
cat("Total genes:  ", nrow(counts_all), "\n")
cat("Total samples:", ncol(counts_all), "\n")
cat("\nFirst 5 genes, first 3 samples:\n")
print(counts_all[1:5, 1:3])


# ============================================================
# STEP 3: DOWNLOAD SAMPLE METADATA FROM GEO
# ============================================================

cat("\n=== DOWNLOADING SAMPLE INFO FROM GEO ===\n")
cat("(This may take 1-2 minutes...)\n")

gse   <- getGEO("GSE66573", GSEMatrix = TRUE)
pheno <- pData(gse[[1]])

cat("Available metadata columns:\n")
print(colnames(pheno))

# Build metadata data frame
metadata_all <- data.frame(
  SampleID         = rownames(pheno),
  Title            = pheno$title,
  stringsAsFactors = FALSE
)

cat("✓ Sample info downloaded\n")
cat("\nAll sample titles:\n")
print(metadata_all$Title)


# ============================================================
# STEP 4: PARSE SAMPLE INFORMATION
# ============================================================

cat("\n=== PARSING SAMPLE INFORMATION ===\n")

# "CTRL" = healthy controls, "RRMS" = relapsing-remitting MS patients
metadata_all$Group <- ifelse(
  grepl("^CTRL", metadata_all$Title),
  "HC",
  "MS"
)

# Verify
print(table(metadata_all$Group))

# Add dataset info
metadata_all$Dataset  <- "GSE66573"
metadata_all$Platform <- "RNA-seq"

cat("✓ Sample information parsed\n")
cat("\nSample breakdown:\n")
print(table(metadata_all$Group))

# Keep only samples that exist in count matrix
metadata_all <- metadata_all[metadata_all$SampleID %in% colnames(counts_all), ]
cat("Samples matching count matrix:", nrow(metadata_all), "\n")

# ============================================================
# STEP 5: QUALITY CONTROL - FILTER LOW-EXPRESSED GENES
# ============================================================
# Small dataset (14 samples) - use lenient threshold
# Keep genes with >= 10 counts in at least 4 samples

cat("\n=== QUALITY CONTROL ===\n")
cat("Genes before filtering:", nrow(counts_all), "\n")

# Subset to our samples first
sample_cols     <- match(metadata_all$SampleID, colnames(counts_all))
counts_filt_pre <- counts_all[, sample_cols]

keep        <- rowSums(counts_filt_pre >= 10) >= 4
counts_filt <- counts_filt_pre[keep, ]

cat("Genes after filtering: ", nrow(counts_filt), "\n")
cat("Genes removed:         ", nrow(counts_filt_pre) - nrow(counts_filt), "\n")


# ============================================================
# STEP 6: CREATE DESeq2 OBJECT
# ============================================================

cat("\n=== CREATING DESeq2 OBJECT ===\n")

# Reorder metadata to match count matrix columns
metadata_all <- metadata_all[match(colnames(counts_filt),
                                   metadata_all$SampleID), ]

# Safety check
stopifnot(all(metadata_all$SampleID == colnames(counts_filt)))

# Make Group a factor with HC as reference
metadata_all$Group <- factor(metadata_all$Group, levels = c("HC", "MS"))

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = metadata_all,
  design    = ~ Group
)

cat("✓ DESeq2 object created\n")
cat("  Genes:  ", nrow(dds), "\n")
cat("  Samples:", ncol(dds), "\n")


# ============================================================
# STEP 7: NORMALISE DATA USING VST
# ============================================================

cat("\n=== NORMALISING DATA (VST) ===\n")

vsd <- vst(dds, blind = TRUE)

cat("✓ Normalisation complete\n")


# ============================================================
# STEP 8: PCA PLOT
# ============================================================

cat("\n=== CREATING PCA PLOT ===\n")

pca_plot <- plotPCA(vsd, intgroup = "Group") +
  ggtitle("GSE66573: MS vs HC") +
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
      "HC" = paste0("Healthy Controls (n=", sum(metadata_all$Group == "HC"), ")"),
      "MS" = paste0("MS Patients (n=",      sum(metadata_all$Group == "MS"), ")")
    )
  )

print(pca_plot)

ggsave(
  filename = "results/figures/PCA_GSE66573.png",
  plot     = pca_plot,
  width    = 7, height = 6, dpi = 300, bg = "white"
)
cat("✓ Saved: results/figures/PCA_GSE66573.png\n")


# ============================================================
# STEP 9: SAVE PROCESSED DATA
# ============================================================

cat("\n=== SAVING PROCESSED DATA ===\n")

saveRDS(counts_filt,  "results/counts_GSE66573_filtered.rds")
saveRDS(metadata_all, "metadata/metadata_GSE66573.rds")
saveRDS(dds,          "results/dds_GSE66573.rds")
saveRDS(vsd,          "results/vsd_GSE66573.rds")

cat("✓ All data saved\n")


# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n")
cat(rep("=", 55), "\n")
cat("             ANALYSIS COMPLETE - GSE66573\n")
cat(rep("=", 55), "\n\n")
cat("Dataset:        GSE66573\n")
cat("Total samples:  ", nrow(metadata_all), "\n")
cat("  MS patients:  ", sum(metadata_all$Group == "MS"), "\n")
cat("  Healthy ctrl: ", sum(metadata_all$Group == "HC"), "\n")
cat("Genes analysed: ", nrow(dds), "\n\n")
cat("Output files:\n")
cat("  results/figures/PCA_GSE66573.png\n")
cat("  results/*.rds\n")
cat("  metadata/*.rds\n")
