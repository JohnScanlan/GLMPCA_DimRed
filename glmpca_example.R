#Sys.setenv(RETICULATE_PYTHON = "C:/Users/crtuser/anaconda/envs/r-reticulate-new/python.exe")

library(reticulate)
use_virtualenv("C:/Users/crtuser/leidenalg-env", required = TRUE)
leidenalg <- reticulate::import("leidenalg")

library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratDisk)
library(glmpca)
library(Matrix)
library(leiden)
library(ggvolc)
library(reticulate)
library(ggplot2)
library(SingleR)
library(dplyr)
library(SingleR)
library(scRNAseq)
library(ExperimentHub)
library(scuttle)

#Patient 1
p1w1.data <- Read10X(data.dir = 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\HumanWeightLossCohort\\WL_PATIENT_10X_Human\\D19-7142\\filtered_feature_bc_matrix')
p1w1 <- CreateSeuratObject(counts = p1w1.data, project = "p1w1")
#p1w3.data <- Read10X(data.dir = 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7143\\filtered_feature_bc_matrix')
#p1w3 <- CreateSeuratObject(counts = p1w3.data, project = "p1w3")#Patient 1
p1w6.data <- Read10X(data.dir = 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\HumanWeightLossCohort\\WL_PATIENT_10X_Human\\D19-7144\\filtered_feature_bc_matrix')
p1w6 <- CreateSeuratObject(counts = p1w6.data, project = "p1w6")#Patient 1

#Patient 2
p2w1.data <- Read10X(data.dir = 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\HumanWeightLossCohort\\WL_PATIENT_10X_Human\\D19-7145\\filtered_feature_bc_matrix')
p2w1 <- CreateSeuratObject(counts = p2w1.data, project = "p2w1")
#p2w3.data <- Read10X(data.dir = 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7146\\filtered_feature_bc_matrix')
#p2w3 <- CreateSeuratObject(counts = p2w3.data, project = "p2w3")#Patient 1
p2w6.data <- Read10X(data.dir = 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\HumanWeightLossCohort\\WL_PATIENT_10X_Human\\D19-7147\\filtered_feature_bc_matrix')
p2w6 <- CreateSeuratObject(counts = p2w6.data, project = "p2w6")

#Patient 3
p3w1.data <- Read10X(data.dir = 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\HumanWeightLossCohort\\WL_PATIENT_10X_Human\\D19-7148\\filtered_feature_bc_matrix')
p3w1 <- CreateSeuratObject(counts = p3w1.data, project = "p3w1")
#p3w3.data <- Read10X(data.dir = 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7149\\filtered_feature_bc_matrix')
#p3w3 <- CreateSeuratObject(counts = p3w3.data, project = "p3w3")
p3w6.data <- Read10X(data.dir = 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\HumanWeightLossCohort\\WL_PATIENT_10X_Human\\D19-7150\\filtered_feature_bc_matrix')
p3w6 <- CreateSeuratObject(counts = p3w6.data, project = "p3w6")

#Patient 1
p4w1.data <- Read10X(data.dir = 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\HumanWeightLossCohort\\WL_PATIENT_10X_Human\\D19-7151\\filtered_feature_bc_matrix')
p4w1 <- CreateSeuratObject(counts = p4w1.data, project = "p4w1")#Patient 1
#p4w3.data <- Read10X(data.dir = 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7152\\filtered_feature_bc_matrix')
#p4w3 <- CreateSeuratObject(counts = p4w3.data, project = "p4w3")#Patient 1
p4w6.data <- Read10X(data.dir = 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\HumanWeightLossCohort\\WL_PATIENT_10X_Human\\D19-7153\\filtered_feature_bc_matrix')
p4w6 <- CreateSeuratObject(counts = p4w6.data, project = "p4w6")

#cells.combined <- readRDS('C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\HumanWeightLossCohort\\WL_PATIENT_10X_Human\\glmpca.rds')

#merge
cells.combined <- merge(p1w1, y=c(p1w6, p2w1, p2w6, p3w1, p3w6, p4w1, p4w6), add.cell.ids = c('p1w1', 'p1w6', 'p2w1', 'p2w6', 'p3w1', 'p3w6', 'p4w1', 'p4w6'), project = 'human weight loss experiment')

# store mitochondrial percentage in object meta data
cells.combined <- PercentageFeatureSet(cells.combined, pattern = "MT-", col.name = "percent.mt")

# store ribosomal percentage in object meta data
cells.combined <- PercentageFeatureSet(cells.combined, pattern = "Rp", col.name = "percent.rp")

# Visualize QC metrics as a violin plot
VlnPlot(cells.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

cells.combined <- subset(cells.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 12.5 & nCount_RNA < 15000)

# run sctransform
cells.combined <- SCTransform(cells.combined, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

# Your existing code for preparing the GLMPCA data
sct_data <- cells.combined@assays$SCT@data
non_zero_rows <- rowSums(sct_data) > 0
sct_data_non_zero <- sct_data[non_zero_rows, ]

# Run GLMPCA on the subsetted data
glmpca_results <- glmpca(Y = sct_data_non_zero, L = 30, minibatch = 'stochastic')

colnames(glmpca_results$factors) <- paste0('glmpca', 1:30)

cells.combined[['glmpca']] <- CreateDimReducObject(embeddings = as.matrix(glmpca_results$factors), loadings = as.matrix(glmpca_results$loadings), assay = DefaultAssay(cells.combined), key = 'glmpca_')

#PCA and UMAP
cells.combined <- RunUMAP(cells.combined, dims = 1:30, verbose = FALSE, reduction = 'glmpca')
DimPlot(cells.combined, group.by = 'orig.ident', reduction = 'glmpca', label = F)

cells.combined <- FindNeighbors(cells.combined, dims = 1:30, verbose = FALSE, reduction = 'glmpca')
cells.combined <- FindClusters(cells.combined, verbose = FALSE, algorithm = 4, method = 'igraph', resolution = 0.3)
DimPlot(cells.combined, reduction = 'umap', label = T)