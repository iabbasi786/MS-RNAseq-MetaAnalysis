# MS RNA-seq Meta-Analysis
**MSc Bioinformatics Research Project**
Iman Abbasi | University of Liverpool

---

## About this project
This is the code for my MSc research project, which looks at gene expression 
differences between Multiple Sclerosis (MS) patients and healthy controls using 
publicly available RNA-seq data from the GEO database.

The initial goal was to run PCA on each dataset to check data quality and see 
whether the datasets could potentially be integrated for a meta-analysis.

---

## Datasets used
All three datasets were downloaded from NCBI GEO:

- **GSE137143** - 137 CD4+ T cell samples (118 MS, 19 HC) sequenced on NovaSeq 6000
- **GSE186895** - 16 whole PBMC samples (8 RIS, 8 HC) sequenced on DNBSEQ-G50
- **GSE66573** - 14 whole blood samples (6 RRMS, 8 HC)

---

## What the scripts do
Each script loads the raw count data, downloads the sample metadata from GEO, 
filters out lowly expressed genes, normalises the data using VST, and produces 
a PCA plot. The scripts are commented throughout to explain each step.

- `GSE137143_PCA_ANALYSIS.R` - includes outlier removal and batch effect investigation
- `GSE186895_PCA_ANALYSIS.R` - RIS vs healthy controls
- `GSE66573_PCA_ANALYSIS.R` - RRMS vs healthy controls

---

## What I found
- GSE137143 showed a clear batch effect after outlier removal
- GSE186895 showed some separation between RIS and HC on PC1
- GSE66573 had outliers and no clear separation between groups
- Batch correction will be needed before any integration can happen

---

## Packages needed
DESeq2, ggplot2, GEOquery - all installation instructions are in the scripts.

---

## Next steps
Batch correction, dataset integration, differential expression analysis 
and pathway enrichment analysis.
