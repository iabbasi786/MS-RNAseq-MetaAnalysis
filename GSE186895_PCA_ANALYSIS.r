# ============================================================
# MS RNA-seq Meta-Analysis
# Script 2: GSE186895 - PCA Analysis
# Dataset: Burlakin et al. 2022 - 16 samples
#          8 Radiologically Isolated Syndrome (RIS) + 8 Healthy Controls
#          Whole PBMCs, DNBSEQ-G50 platform
# Note: RIS = pre-clinical MS (first demyelinating lesion,
#       no clinical symptoms yet). Used as proxy for early MS.
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
# Your working directory should be the GSE186895 data folder
# Change this to wherever you have stored the GSE186895 files

cat("Working directory:", getwd(), "\n")

# Create output folders
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("metadata",        recursive = TRUE, showWarnings = FALSE)


# ============================================================
# STEP 2: LOAD THE COUNT MATRIX
# ============================================================
# Small dataset: 16 samples, whole PBMCs
# Rows = genes, Columns = samples

counts_all <- read.table(
  file         = "GSE186895_raw_counts_GRCh38.p13_NCBI.tsv",
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

gse   <- getGEO("GSE186895", GSEMatrix = TRUE)
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
# Check what the titles look like before parsing
# then extract Group accordingly

cat("\n=== PARSING SAMPLE INFORMATION ===\n")

# Look at titles to understand naming convention
cat("Sample titles:\n")
print(pheno$title)

# Extract Group
# Titles starting with "C_" = healthy controls
# Titles starting with "RIS_" = Radiologically Isolated Syndrome patients
metadata_all$Group <- ifelse(
  grepl("^C_", metadata_all$Title),
  "HC",
  "RIS"
)

# Verify it worked
print(table(metadata_all$Group))

# Add dataset info
metadata_all$Dataset  <- "GSE186895"
metadata_all$Platform <- "DNBSEQ-G50"

cat("✓ Sample information parsed\n")
cat("\nSample breakdown:\n")
print(table(metadata_all$Group))

# Keep only samples that exist in count matrix
metadata_all <- metadata_all[metadata_all$SampleID %in% colnames(counts_all), ]
cat("Samples matching count matrix:", nrow(metadata_all), "\n")


# ============================================================
# STEP 5: QUALITY CONTROL - FILTER LOW-EXPRESSED GENES
# ============================================================
# With only 16 samples, we use a more lenient threshold
# Rule: keep genes with >= 10 counts in at least 5 samples
# (using 10 samples would be too strict for a 16-sample dataset)

cat("\n=== QUALITY CONTROL ===\n")
cat("Genes before filtering:", nrow(counts_all), "\n")

# Subset count matrix to match metadata samples
sample_cols  <- match(metadata_all$SampleID, colnames(counts_all))
counts_filt_pre <- counts_all[, sample_cols]

keep         <- rowSums(counts_filt_pre >= 10) >= 5
counts_filt  <- counts_filt_pre[keep, ]

cat("Genes after filtering: ", nrow(counts_filt), "\n")
cat("Genes removed:         ", nrow(counts_filt_pre) - nrow(counts_filt), "\n")


# ============================================================
# STEP 6: CREATE DESeq2 OBJECT
# ============================================================

cat("\n=== CREATING DESeq2 OBJECT ===\n")

# Reorder metadata to match filtered count matrix columns
metadata_all <- metadata_all[match(colnames(counts_filt),
                                   metadata_all$SampleID), ]

# Safety check
stopifnot(all(metadata_all$SampleID == colnames(counts_filt)))

# Make Group a factor with HC as reference
metadata_all$Group <- factor(metadata_all$Group, levels = c("HC", "RIS"))

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
  ggtitle("GSE186895: RIS vs HC (Whole PBMCs)") +
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
    values = c("HC" = "#3498DB", "RIS" = "#E67E22"),
    name   = "Group",
    labels = c(
      "HC"  = paste0("Healthy Controls (n=", sum(metadata_all$Group == "HC"),  ")"),
      "RIS" = paste0("RIS Patients (n=",     sum(metadata_all$Group == "RIS"), ")")
    )
  )

print(pca_plot)

ggsave(
  filename = "results/figures/PCA_GSE186895.png",
  plot     = pca_plot,
  width    = 7, height = 6, dpi = 300, bg = "white"
)
cat("✓ Saved: results/figures/PCA_GSE186895.png\n")


# ============================================================
# STEP 9: SAVE PROCESSED DATA
# ============================================================

cat("\n=== SAVING PROCESSED DATA ===\n")

saveRDS(counts_filt,  "results/counts_GSE186895_filtered.rds")
saveRDS(metadata_all, "metadata/metadata_GSE186895.rds")
saveRDS(dds,          "results/dds_GSE186895.rds")
saveRDS(vsd,          "results/vsd_GSE186895.rds")

cat("✓ All data saved\n")


# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n")
cat(rep("=", 55), "\n")
cat("             ANALYSIS COMPLETE - GSE186895\n")
cat(rep("=", 55), "\n\n")
cat("Dataset:        GSE186895 (Whole PBMCs)\n")
cat("Total samples:  ", nrow(metadata_all), "\n")
cat("  RIS patients: ", sum(metadata_all$Group == "RIS"), "\n")
cat("  Healthy ctrl: ", sum(metadata_all$Group == "HC"),  "\n")
cat("Genes analysed: ", nrow(dds), "\n\n")
cat("Output files:\n")
cat("  results/figures/PCA_GSE186895.png\n")
cat("  results/*.rds\n")
cat("  metadata/*.rds\n\n")
cat("Note: RIS (Radiologically Isolated Syndrome) = pre-clinical MS.\n")
cat("Interpret separation from HC with this caveat in mind.\n")
