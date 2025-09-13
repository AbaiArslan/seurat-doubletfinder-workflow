# filter out doublets: DoubletFinder

### load libraries----
library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)

### Unzip .gz file----
tar("data/10k_PBMC_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix.tar.gz", exdir = "PBMC_matrix")

### create counts matrix----
cts <- ReadMtx(mtx = 'data/raw_feature_bc_matrix/matrix.mtx.gz',
        features ='data/raw_feature_bc_matrix/features.tsv.gz',
        cells ='data/raw_feature_bc_matrix/barcodes.tsv.gz')


### Create Seurat Object----
pbmc_seurat <- CreateSeuratObject(counts = cts)

### QC and Filtering----
#explore QC
view(pbmc_seurat@meta.data)
#mitondrial percentage
pbmc_seurat[["percent_mt"]] <- PercentageFeatureSet(pbmc_seurat, pattern = "^MT-")
view(pbmc_seurat@meta.data) #viewing data after adding MT Column
VlnPlot(pbmc_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(pbmc_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# filtering
pbmc_seurat_filtering <- subset(pbmc_seurat, subset = nCount_RNA > 800 &
         nFeature_RNA > 500 &
         percent_mt < 10)
view (pbmc_seurat_filtering)
pbmc_fs <- pbmc_seurat_filtering

# pre-process standard workflow
pbmc_fs <- NormalizeData(object = pbmc_fs)
pbmc_fs <- FindVariableFeatures(pbmc_fs)
pbmc_fs <- ScaleData(object = pbmc_fs)
pbmc_fs <- RunPCA(object = pbmc_fs)
ElbowPlot(pbmc_fs)
pbmc_fs <- FindNeighbors(object = pbmc_fs, dims = 1:20)
pbmc_fs <- FindClusters(object = pbmc_fs)
pbmc_fs <- RunUMAP(object = pbmc_fs, dims = 1:20)


### pK identification (no_ground_truth)----
sweep_pbmc_fs <- paramSweep(pbmc_fs, PCs = 1:20, sct = FALSE)
sum_pbmc_fs <- summarizeSweep(sweep_pbmc_fs, GT = F)
bcmvn_pbmc <- find.pK(sum_pbmc_fs)

ggplot(bcmvn_pbmc, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_pbmc %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

### Hymotypic Doublets----

annotations <- pbmc_fs@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*nrow(pbmc_fs@meta.data))
nExp_poi_adj <- round(nExp_poi*(1-homotypic.prop))


### Doublet finder----
library(DoubletFinder)
pbmc_fs_dup <- DoubletFinder::doubletFinder(pbmc_fs,
                         PCs = 1:20,
                         pN = 0.25,
                         pK = pK,
                         nExp = nExp_poi,
                         reuse.pANN = F, sct = F)