# GLMPCA_DimRed
Dimensionality Reduction Using GLMPCA

# Human Weight Loss Cohort Analysis

This repository contains code and instructions for analyzing single-cell RNA-seq data from a human weight loss cohort. The analysis leverages the Seurat R package for preprocessing, quality control, clustering, and dimensionality reduction, as well as custom scripts for GLMPCA-based analysis.

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Data Preparation](#data-preparation)
- [Pipeline Steps](#pipeline-steps)
- [Usage](#usage)
- [License](#license)

---

## Overview

This project processes and analyzes scRNA-seq data from multiple time points (Week 1 and Week 6) for patients undergoing weight loss interventions. Key functionalities include:

- Loading and preprocessing raw 10X Genomics data.
- Performing quality control (QC) to filter cells based on mitochondrial and ribosomal content.
- Data normalization using SCTransform.
- Dimensionality reduction with GLMPCA and UMAP.
- Clustering and differential expression analysis.

---

## Requirements

- R version 4.0 or higher
- The following R packages:
  - `reticulate`
  - `Seurat`
  - `dplyr`
  - `ggplot2`
  - `SeuratDisk`
  - `glmpca`
  - `Matrix`
  - `leiden`
  - `ggvolc`
  - `SingleR`
  - `scRNAseq`
  - `ExperimentHub`
  - `scuttle`
- Python environment with `leidenalg` installed.

---

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/human-weight-loss-cohort.git
   cd human-weight-loss-cohort
   ```

2. Install R packages:
   ```R
   install.packages(c("reticulate", "Seurat", "dplyr", "ggplot2", "SeuratDisk", "glmpca", "Matrix", "leiden"))
   BiocManager::install(c("SingleR", "scRNAseq", "ExperimentHub", "scuttle"))
   ```

3. Set up the Python environment for `leidenalg`:
   ```R
   library(reticulate)
   use_virtualenv("C:/Users/crtuser/leidenalg-env", required = TRUE)
   ```

---

## Data Preparation

1. Prepare your 10X Genomics data directories for each patient and time point.
2. Ensure the paths to these directories are correctly specified in the script.

Example:
```R
p1w1.data <- Read10X(data.dir = 'C:/path/to/patient1_week1/')
```

---

## Pipeline Steps

1. **Data Loading:**
   - Load raw 10X Genomics data using `Read10X`.
   - Create Seurat objects for each sample.

2. **Merging Data:**
   - Combine all samples into one Seurat object using `merge()`.

3. **Quality Control:**
   - Visualize QC metrics with `VlnPlot`.
   - Filter cells based on feature count, mitochondrial percentage, and RNA count thresholds.

4. **Normalization:**
   - Use SCTransform to normalize and scale the data.

5. **Dimensionality Reduction:**
   - Run GLMPCA on normalized data.
   - Generate UMAP visualizations.

6. **Clustering:**
   - Perform clustering using the Leiden algorithm.

7. **Visualization:**
   - Generate plots such as `DimPlot` for cluster visualization.

---

## Usage

Run the script sequentially to reproduce the analysis. Adjust file paths, thresholds, and parameters as needed.

### Key Commands

1. **Run SCTransform and GLMPCA:**
   ```R
   cells.combined <- SCTransform(cells.combined, method = "glmGamPoi", vars.to.regress = "percent.mt")
   
   sct_data <- cells.combined@assays$SCT@data
   non_zero_rows <- rowSums(sct_data) > 0
   sct_data_non_zero <- sct_data[non_zero_rows, ]

   glmpca_results <- glmpca(Y = sct_data_non_zero, L = 30, minibatch = 'stochastic')
   cells.combined[['glmpca']] <- CreateDimReducObject(embeddings = as.matrix(glmpca_results$factors), loadings = as.matrix(glmpca_results$loadings))
   ```

2. **Run UMAP and Clustering:**
   ```R
   cells.combined <- RunUMAP(cells.combined, dims = 1:30, reduction = 'glmpca')
   cells.combined <- FindNeighbors(cells.combined, dims = 1:30, reduction = 'glmpca')
   cells.combined <- FindClusters(cells.combined, resolution = 0.3)
   ```

3. **Visualization:**
   ```R
   DimPlot(cells.combined, group.by = 'orig.ident', label = TRUE)
   ```

---

## License

This project is licensed under the MIT License. See the LICENSE file for details.
