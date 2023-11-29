library(imcRtools)
library(cytomapper)
library(openxlsx)
library(stringr)
library(dittoSeq)
library(RColorBrewer)

#7- Batch correction

#7.1- fastMNN correction with Batchelor
# 7.1.1- Perform sample correction
library(batchelor)
set.seed(220228)
out <- fastMNN(spe, batch = spe$patient_id,
               auto.merge = TRUE,
               subset.row = rowData(spe)$use_channel,
               assay.type = "exprs")

# 7.1.2- Transfer the correction results to the main spe object
reducedDim(spe, "fastMNN") <- reducedDim(out, "corrected")

# 7.1.3- Quality control of correction results
#The fastMNN function further returns outputs that can be used to assess the quality of the batch correction. 
#The metadata(out)$merge.info entry collects diagnostics for each individual merging step. 
#Here, the batch.size and lost.var entries are important. 
#The batch.size entry reports the relative magnitude of the batch effect 
#and the lost.var entry represents the percentage of lost variance per merging step. 
#A large batch.size and low lost.var indicate sufficient batch correction.
merge_info <- metadata(out)$merge.info 

DataFrame(left = merge_info$left,
          right = merge_info$right,
          batch.size = merge_info$batch.size,
          max_lost_var = rowMax(merge_info$lost.var))


# 7.1.4- Visualization
#The simplest option to check if the sample effects were corrected is by using non-linear dimensionality reduction techniques 
#and observe mixing of cells across samples. 
#We will recompute the UMAP embedding using the corrected low-dimensional coordinates for each cell.
library(scater)
set.seed(220228)
spe <- runUMAP(spe, dimred= "fastMNN", name = "UMAP_mnnCorrected") 

#Visualize the corrected UMAP while overlaying patient IDs.
library(cowplot)
library(dittoSeq)
library(viridis)

# visualize patient id 
p1 <- dittoDimPlot(spe, var = "patient_id", 
                   reduction.use = "UMAP", size = 0.2) + 
  scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
  ggtitle("Patient ID on UMAP before correction")
p2 <- dittoDimPlot(spe, var = "patient_id", 
                   reduction.use = "UMAP_mnnCorrected", size = 0.2) + 
  scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
  ggtitle("Patient ID on UMAP after correction")

plot_grid(p1, p2)

#visualize the expression of selected markers across all cells before batch correction.
markers <- c("Ecad", "CD45RO", "CD20", "CD3", "FOXP3", "CD206", "MPO", "SMA", "Ki67")

# Before correction
plot_list <- multi_dittoDimPlot(spe, var = markers, reduction.use = "UMAP", 
                                assay = "exprs", size = 0.2, list.out = TRUE) 
plot_list <- lapply(plot_list, function(x) x + scale_color_viridis())
plot_grid(plotlist = plot_list) 

##visualize the expression of selected markers across all cells after batch correction.
# After correction
plot_list <- multi_dittoDimPlot(spe, var = markers, reduction.use = "UMAP_mnnCorrected", 
                                assay = "exprs", size = 0.2, list.out = TRUE) 
plot_list <- lapply(plot_list, function(x) x + scale_color_viridis())
plot_grid(plotlist = plot_list) 


#7.2- Harmony
# First create the expression matrix and call the HarmonyMatrix function to perform the correction.
library(harmony)
mat <- t(assay(spe, "exprs")[rowData(spe)$use_channel,])
harmony_emb <- HarmonyMatrix(mat, spe$patient_id, do_pca = TRUE)
reducedDim(spe, "harmony") <- harmony_emb

#7.2.1- Visualization
set.seed(220228)
spe <- runUMAP(spe, dimred = "harmony", name = "UMAP_harmony") 

# visualize patient id 
p1 <- dittoDimPlot(spe, var = "patient_id", 
                   reduction.use = "UMAP", size = 0.2) + 
  scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
  ggtitle("Patient ID on UMAP before correction")
p2 <- dittoDimPlot(spe, var = "patient_id", 
                   reduction.use = "UMAP_harmony", size = 0.2) + 
  scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
  ggtitle("Patient ID on UMAP after correction")

plot_grid(p1, p2)

# Before correction
plot_list <- multi_dittoDimPlot(spe, var = markers, reduction.use = "UMAP", 
                                assay = "exprs", size = 0.2, list.out = TRUE) 
plot_list <- lapply(plot_list, function(x) x + scale_color_viridis())
plot_grid(plotlist = plot_list) 

# After correction
plot_list <- multi_dittoDimPlot(spe, var = markers, reduction.use = "UMAP_harmony", 
                                assay = "exprs", size = 0.2, list.out = TRUE) 
plot_list <- lapply(plot_list, function(x) x + scale_color_viridis())
plot_grid(plotlist = plot_list) 

#7.3- Seurat correction
#To use Seurat, we will first create a Seurat object from the SpatialExperiment object and add relevant metadata. 
#The object also needs to be split by patient prior to integration.
library(Seurat)
library(SeuratObject)
seurat_obj <- as.Seurat(spe, counts = "counts", data = "exprs")
seurat_obj <- AddMetaData(seurat_obj, as.data.frame(colData(spe)))
seurat.list <- SplitObject(seurat_obj, split.by = "patient_id")

#To avoid long run times, we will use an approach that relies on reciprocal PCA 
#instead of canonical correlation analysis for dimensionality reduction and initial alignment.
#We will first define the features used for integration and perform PCA on cells of each patient individually. 
features <- rownames(spe)[rowData(spe)$use_channel]
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  return(x)
})

#The FindIntegrationAnchors function detects MNNs between cells of different patients 
#and the IntegrateData function corrects the expression values of cells. 
#We slightly increase the number of neighbors to be considered for MNN detection (the k.anchor parameter). 
#This increases the integration strength.
anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                  anchor.features = features,
                                  reduction = "rpca", 
                                  k.anchor = 20)

combined <- IntegrateData(anchorset = anchors)

#We now select the integrated assay and perform PCA dimensionality reduction. 
#The cell coordinates in PCA reduced space can then be transferred to the original SpatialExperiment object.
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
reducedDim(spe, "seurat") <- combined@reductions$pca@cell.embeddings

# 7.3.1- Visualization
# We recompute the UMAP embeddings based on Seurat integrated results and visualize the embedding.
set.seed(220228)
spe <- runUMAP(spe, dimred = "seurat", name = "UMAP_seurat") 

# visualize patient id 
p1 <- dittoDimPlot(spe, var = "patient_id", 
                   reduction.use = "UMAP", size = 0.2) + 
  scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
  ggtitle("Patient ID on UMAP before correction")
p2 <- dittoDimPlot(spe, var = "patient_id", 
                   reduction.use = "UMAP_seurat", size = 0.2) + 
  scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
  ggtitle("Patient ID on UMAP after correction")

plot_grid(p1, p2)

# Before correction
plot_list <- multi_dittoDimPlot(spe, var = markers, reduction.use = "UMAP", 
                                assay = "exprs", size = 0.2, list.out = TRUE) 
plot_list <- lapply(plot_list, function(x) x + scale_color_viridis())
plot_grid(plotlist = plot_list) 

# After correction
plot_list <- multi_dittoDimPlot(spe, var = markers, reduction.use = "UMAP_seurat", 
                                assay = "exprs", size = 0.2, list.out = TRUE) 
plot_list <- lapply(plot_list, function(x) x + scale_color_viridis())
plot_grid(plotlist = plot_list) 

#The modified SpatialExperiment object is saved for further downstream analysis.
saveRDS(spe, "spe.rds")
spe <- readRDS("spe.rds")
