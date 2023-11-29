library(imcRtools)
library(cytomapper)
library(openxlsx)
library(stringr)
library(dittoSeq)
library(RColorBrewer)
library(Rphenograph)
library(igraph)
library(dittoSeq)
library(viridis)
library(bluster)
library(BiocParallel)
library(ggplot2)
library(scran)
library(CATALYST)
library(kohonen)
library(ConsensusClusterPlus)
library(patchwork)
library(pheatmap)
library(gridExtra)
library(SingleCellExperiment)
library(tidyverse)
library(ggridges)

#1- Load data
#First, we will read in the previously generated SpatialExperiment object.

spe <- readRDS("spe.rds")

#For visualization purposes, we will define markers that were used for cell type classification 
#and markers that can indicate a specific cell state (e.g. Ki67).

# Define cell_type_markers 
type_markers <- c("Ecad", "CD45RO", "CD20", "CD3", "FOXP3", "CD206", "MPO", 
                  "SMA", "CD8a", "CD4", "HLADR", "CD15", "CD38", "PDGFRb")

# Define cell_state_markers 
state_markers <- c("CarbonicAnhydrase", "Ki67", "PD1", "GrzB", "PDL1", 
                   "ICOS", "TCF7", "VISTA")

# Add to spe
rowData(spe)$marker_class <- ifelse(rownames(spe) %in% type_markers, "type",
                                    ifelse(rownames(spe) %in% state_markers, "state", 
                                           "other"))

#2- Cell-type level

#In the first section of this chapter, the grouping-level for the visualization approaches 
#will be the cell type classification from the previous script. 
#Other grouping levels (e.g. cluster assignments from Section 9.2) are possible and 
#the user should adjust depending on the chosen analysis workflow.

#2.1- Dimensionality reduction visualization

#As seen before, we can visualize single-cells in low-dimensional space. 
#Often, non-linear methods for dimensionality reduction such as tSNE and UMAP are sued. 
#They aim to preserve the distances between each cell and its neighbors in the high-dimensional space.
#Here, we will use dittoDimPlot from the DittoSeq package and plotReducedDim from the scater package 
#to visualize the fastMNN-corrected UMAP colored by cell type and expression, respectively.
#Both functions are highly flexible and return ggplot objects which can be further modified.

library(dittoSeq)
library(scater)
library(patchwork)
library(cowplot)
library(viridis)


## UMAP colored by cell type and expression - dittoDimPlot
p1 <- dittoDimPlot(spe, var = "celltype",
                   reduction.use = "UMAP_mnnCorrected", size = 0.2,
                   do.label = TRUE) +
  scale_color_manual(values = metadata(spe)$color_vectors$celltype) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell types on UMAP, integrated cells")


p2 <- dittoDimPlot(spe, var = "Ecad", assay = "exprs",
                   reduction.use = "UMAP_mnnCorrected", size = 0.2, 
                   colors = viridis(100), do.label = TRUE) +
  scale_color_viridis()

p1 + p2

# UMAP colored by expression for all markers - plotReducedDim
plot_list  <- lapply(rownames(spe)[rowData(spe)$marker_class == "type"], function(x){
  p <- plotReducedDim(spe, dimred = "UMAP_mnnCorrected",
                      colour_by = x,
                      by_exprs_values = "exprs",
                      point_size = 0.2)
  return(p)
})

plot_grid(plotlist = plot_list)

#2.2- Heatmap visualisation

#ext, it is often useful to visualize single-cell expression per cell type in form of a heatmap. 
#For this, we will use the dittoHeatmap function from the DittoSeq package.
#We sub-sample the dataset to 4000 cells for ease of visualization and 
#overlay the cancer type and patient ID from which the cells were extracted.

set.seed(220818)
cur_cells <- sample(seq_len(ncol(spe)), 4000)

#Heatmap visualization - DittoHeatmap
dittoHeatmap(spe[,cur_cells], genes = rownames(spe)[rowData(spe)$marker_class == "type"],
             assay = "exprs", order.by = c("celltype"),
             cluster_cols = FALSE, scale = "none",
             heatmap.colors = viridis(100), annot.by = c("celltype","indication","patient_id"),
             annotation_colors = list(indication = metadata(spe)$color_vectors$indication,
                                      patient_id = metadata(spe)$color_vectors$patient_id,
                                      celltype = metadata(spe)$color_vectors$celltype)
)
#annotation_colors = list(indication = metadata(spe)$color_vectors$indication,patient_id = metadata(spe)$color_vectors$patient_id)

#Similarly, we can visualize the mean marker expression per cell type for all cells using aggregateAcrossCells from scuttle 
#and then use dittoHeatmap. We will annotate the heatmap with the number of cells per cell type.

library(scuttle)
## by cell type
celltype_mean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),  
                                      ids = spe$celltype, 
                                      statistics = "mean",
                                      use.assay.type = "exprs", 
                                      subset.row = rownames(spe)[rowData(spe)$marker_class == "type"]
)

# No scaling
dittoHeatmap(celltype_mean,
             assay = "exprs", cluster_cols = TRUE, 
             scale = "none",
             heatmap.colors = viridis(100),
             annot.by = c("celltype","ncells"),
             annotation_colors = list(celltype = metadata(spe)$color_vectors$celltype,
                                      ncells = plasma(100)))

# Min-max expression scaling
dittoHeatmap(celltype_mean,
             assay = "exprs", cluster_cols = TRUE, 
             scaled.to.max = TRUE,
             heatmap.colors.max.scaled = inferno(100),
             annot.by = c("celltype","ncells"),
             annotation_colors = list(celltype = metadata(spe)$color_vectors$celltype,
                                      ncells = plasma(100)))

#2.3-  Violin plot visualization

#The plotExpression function from the scater package allows to plot the distribution of expression values 
#across cell types for a chosen set of proteins. The output is a flexible ggplot object.

#Violin Plot - plotExpression
plotExpression(spe[,cur_cells], 
               features = rownames(spe)[rowData(spe)$marker_class == "type"],
               x = "celltype", exprs_values = "exprs", 
               colour_by = "celltype") +
  theme(axis.text.x =  element_text(angle = 90))+
  scale_color_manual(values = metadata(spe)$color_vectors$celltype)

#2.4- Scatter plot visualization

#Moreover, a protein expression based scatter plot can be generated with dittoScatterPlot (returns a ggplot object). 
#We overlay the plot with the cell type information.

#Scatter plot
dittoScatterPlot(spe, 
                 x.var = "CD3", y.var="CD20", 
                 assay.x = "exprs", assay.y = "exprs", 
                 color.var = "celltype") +
  scale_color_manual(values = metadata(spe)$color_vectors$celltype) +
  ggtitle("Scatterplot for CD3/CD20 labelled by celltype")


#2.5- Barplot visualization

# by sample_id
dittoBarPlot(spe, var = "pg_clusters", group.by = "sample_id") +
  scale_fill_manual(values = c(dittoColors(1)[1:length(unique(spe$pg_clusters))]))

# by patient_id - percentage
# by sample_id
# by patient_id - percentage
dittoBarPlot(spe, var = "celltype", group.by = "patient_id") +
  scale_fill_manual(values = metadata(spe)$color_vectors$celltype)
# by patient_id - count
dittoBarPlot(spe, scale = "count", var = "celltype", group.by = "patient_id") +
  scale_fill_manual(values = metadata(spe)$color_vectors$celltype)

# by disease group - percentage
dittoBarPlot(spe, var = "celltype", group.by = "indication") +
  scale_fill_manual(values = metadata(spe)$color_vectors$celltype)

#2.6- CATALYST-based visualization
#In the following, we highlight some useful visualization functions from the CATALYST package.
#To this end, we will first convert the SpatialExperiment object into a CATALYST-compatible format.

library(CATALYST)
# save spe in CATALYST-compatible object with renamed colData entries and 
# new metadata information
spe_cat <- spe 

spe_cat$sample_id <- factor(spe$sample_id)
spe_cat$condition <- factor(spe$indication)
spe_cat$cluster_id <- factor(spe$celltype)

#add celltype information to metadata
metadata(spe_cat)$cluster_codes <- data.frame(celltype = factor(spe_cat$cluster_id))

#2.6.1- Pseudobulk-level MDS plot
#Pseudobulk-level multi-dimensional scaling (MDS) plots can be rendered with the exported pbMDS function.
#Here, we will use pbMDS to highlight expression similarities between cell types and subsequently for each celltype-sample-combination

# MDS pseudobulk by cell type
pbMDS(spe_cat, by = "cluster_id", 
      features = rownames(spe_cat)[rowData(spe_cat)$marker_class == "type"], 
      label_by = "cluster_id", k = "celltype") +
  scale_color_manual(values = metadata(spe_cat)$color_vectors$celltype)

# MDS pseudobulk by cell type and sample_id
pbMDS(spe_cat, by = "both", 
      features = rownames(spe_cat)[rowData(spe_cat)$marker_class == "type"], 
      k = "celltype", shape_by = "condition", 
      size_by = TRUE) +
  scale_color_manual(values = metadata(spe_cat)$color_vectors$celltype)

#2.6.2- Reduced dimension plot on CLR of proportions
#The clrDR function produces dimensionality reduction plots on centered log-ratios (CLR) of sample/cell type proportions across cell type/samples.
#As with pbMDS, the output plots aim to illustrate the degree of similarity between cell types based on sample proportions.

# CLR on cluster proportions across samples
clrDR(spe_cat, dr = "PCA", 
      by = "cluster_id", k = "celltype", 
      label_by = "cluster_id", arrow_col = "patient_id", 
      point_pal = metadata(spe_cat)$color_vectors$celltype) +
  scale_color_manual(values = metadata(spe_cat)$color_vectors$patient_id)

#2.6.3- Pseudobulk expression boxplot
#The plotPbExprs generates combined box- and jitter-plots of aggregated marker expression per cell type. 
#Here, we further split the data by cancer type.
plotPbExprs(spe_cat, k = "celltype", 
            facet_by = "cluster_id", ncol = 4, 
            features = rownames(spe_cat)[rowData(spe_cat)$marker_class == "type"]) +
  scale_color_manual(values = metadata(spe_cat)$color_vectors$indication)


#3- Sample-level

#3.1- Dimensionality reduction visualization
#Visualization of low-dimensional embeddings, here comparing non-corrected and fastMNN-corrected UMAPs, 
#and coloring it by sample-levels is often used for “batch effect” assessment as mentioned in Section 7.4.
#We will again use dittoDimPlot.

## UMAP colored by cell type and expression - dittoDimPlot
p1 <- dittoDimPlot(spe, var = "sample_id",
                   reduction.use = "UMAP", size = 0.2, 
                   colors = viridis(100), do.label = FALSE) +
  scale_color_manual(values = metadata(spe)$color_vectors$sample_id) +
  theme(legend.title = element_blank()) +
  ggtitle("Sample ID")

p2 <- dittoDimPlot(spe, var = "sample_id",
                   reduction.use = "UMAP_mnnCorrected", size = 0.2, 
                   colors = viridis(100), do.label = FALSE) +
  scale_color_manual(values = metadata(spe)$color_vectors$sample_id) +
  theme(legend.title = element_blank()) +
  ggtitle("Sample ID")

p3 <- dittoDimPlot(spe, var = "patient_id", 
                   reduction.use = "UMAP", size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
  theme(legend.title = element_blank()) +
  ggtitle("Patient ID")

p4 <- dittoDimPlot(spe, var = "patient_id", 
                   reduction.use = "UMAP_mnnCorrected", size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
  theme(legend.title = element_blank()) +
  ggtitle("Patient ID")

(p1 + p2) / (p3 + p4)

#3.2- Heatmap visualization
#It can be beneficial to use a heatmap to visualize single-cell expression per sample and patient. 
#Such a plot, which we will create using dittoHeatmap, can highlight biological differences across samples/patients.
#Heatmap visualization - DittoHeatmap
dittoHeatmap(spe[,cur_cells], genes = rownames(spe)[rowData(spe)$marker_class == "type"],
             assay = "exprs", order.by = c("patient_id","sample_id"),
             cluster_cols = FALSE, scale = "none",
             heatmap.colors = viridis(100), 
             annot.by = c("celltype","indication","patient_id","sample_id"),
             annotation_colors = list(celltype = metadata(spe)$color_vectors$celltype,
                                      indication = metadata(spe)$color_vectors$indication,
                                      patient_id = metadata(spe)$color_vectors$patient_id,
                                      sample_id = metadata(spe)$color_vectors$sample_id))

#aggregated mean marker expression per sample/patient allow identification of samples/patients with outlying expression patterns.
#Here, we will focus on the patient level and use aggregateAcrossCells and dittoHeatmap. 
#The heatmap will be annotated with the number of cells per patient and cancer type and displayed using two scaling options.

#by patient_id
patient_mean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),  
                                     ids = spe$patient_id, 
                                     statistics = "mean",
                                     use.assay.type = "exprs", 
                                     subset.row = rownames(spe)[rowData(spe)$marker_class == "type"]
)

# No scaling
dittoHeatmap(patient_mean,
             assay = "exprs", cluster_cols = TRUE, 
             scale = "none",
             heatmap.colors = viridis(100),
             annot.by = c("patient_id","indication","ncells"),
             annotation_colors = list(patient_id = metadata(spe)$color_vectors$patient_id,
                                      indication = metadata(spe)$color_vectors$indication,
                                      ncells = plasma(100)))

# Min-max expression scaling
dittoHeatmap(patient_mean,
             assay = "exprs", cluster_cols = TRUE, 
             scaled.to.max =  TRUE,
             heatmap.colors.max.scaled = inferno(100),
             annot.by = c("patient_id","indication","ncells"),
             annotation_colors = list(patient_id = metadata(spe)$color_vectors$patient_id,
                                      indication = metadata(spe)$color_vectors$indication,
                                      ncells = plasma(100)))

#3.3- Barplot visualization

#Complementary to displaying cell type frequencies per sample/patient, 
#we can use dittoBarPlot to display sample/patient frequencies per cell type.
dittoBarPlot(spe, var = "patient_id", group.by = "pg_clusters") +
  scale_fill_manual(values = metadata(spe)$color_vectors$patient_id)

dittoBarPlot(spe, var = "sample_id", group.by = "pg_clusters") +
  scale_fill_manual(values = metadata(spe)$color_vectors$sample_id)


#3.4- CATALYST-based visualization

#3.4.1 Pseudobulk-level MDS plot
#Expression-based pseudobulks for each sample can be compared with the pbMDS function.


# MDS pseudobulk by sample_id 
pbMDS(spe_cat, by = "sample_id", 
      color_by = "sample_id", 
      features = rownames(spe_cat)[rowData(spe_cat)$marker_class == "type"]) +
  scale_color_manual(values = metadata(spe_cat)$color_vectors$sample_id)


#3.4.2- Reduced dimension plot on CLR of proportions

#The clrDR function can also be used to analyze similarity of samples based on cell type proportions.
# CLR on sample proportions across clusters
clrDR(spe_cat, dr = "PCA", 
      by = "sample_id", point_col = "sample_id",
      k = "celltype", point_pal = metadata(spe_cat)$color_vectors$sample_id) +
  scale_color_manual(values = metadata(spe_cat)$color_vectors$celltype)

#4- Publication-ready ComplexHeatmap

#For this example, we will concatenate heatmaps and annotations horizontally into one rich heatmap list. The grouping-level for the visualization will again be the cell type information from Section 9.3
#Initially, we will create two separate Heatmap objects for cell type and state markers.
#Then, metadata information, including the cancer type proportion and number of cells/patients per cell type, will be extracted into HeatmapAnnotation objects.
#Notably, we will add spatial features per cell type, here the number of neighbors extracted from colPair(spe) and cell area, in another HeatmapAnnotation object.
#Ultimately, all objects are combined in a HeatmapList and visualized.
### 1. Heatmap bodies ###

library(ComplexHeatmap)
library(circlize)
library(tidyverse)
set.seed(22)

# Heatmap body color 
col_exprs <- colorRamp2(c(0,1,2,3,4), 
                        c("#440154FF","#3B518BFF","#20938CFF",
                          "#6ACD5AFF","#FDE725FF"))

# Create Heatmap objects
# By celltype markers
celltype_mean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),  
                                      ids = spe$celltype, 
                                      statistics = "mean",
                                      use.assay.type = "exprs", 
                                      subset.row = rownames(spe)[rowData(spe)$marker_class == "type"])

h_type <- Heatmap(t(assay(celltype_mean, "exprs")),
                  column_title = "type_markers",
                  col = col_exprs,
                  name= "mean exprs",
                  show_row_names = TRUE, 
                  show_column_names = TRUE)

# By cellstate markers
cellstate_mean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),  
                                       ids = spe$celltype, 
                                       statistics = "mean",
                                       use.assay.type = "exprs", 
                                       subset.row = rownames(spe)[rowData(spe)$marker_class == "state"])

h_state <- Heatmap(t(assay(cellstate_mean, "exprs")),
                   column_title = "state_markers",
                   col = col_exprs,
                   name= "mean exprs",
                   show_row_names = TRUE,
                   show_column_names = TRUE)


### 2. Heatmap annotation ###

### 2.1  Metadata features

anno <- colData(celltype_mean) %>% as.data.frame %>% select(celltype, ncells)

# Proportion of indication per celltype
indication <- colData(spe) %>% 
  as.data.frame() %>% 
  select(celltype, indication) %>% 
  group_by(celltype) %>% 
  table() %>% 
  as.data.frame()

indication <- indication %>% 
  group_by(celltype) %>% 
  mutate(fra = Freq/sum(Freq)) 

indication <- indication %>% 
  select(-Freq) %>% 
  pivot_wider(id_cols = celltype, 
              names_from = indication, 
              values_from = fra) %>% 
  column_to_rownames("celltype")

# Number of contributing patients per celltype
cluster_PID <- colData(spe) %>% 
  as.data.frame() %>% 
  select(celltype, patient_id) %>% 
  group_by(celltype) %>% table() %>% 
  as.data.frame()

n_PID <- cluster_PID %>% 
  filter(Freq>0) %>% 
  group_by(celltype) %>% 
  count(name = "n_PID") %>% 
  column_to_rownames("celltype")

# Create HeatmapAnnotation objects
ha_anno <- HeatmapAnnotation(celltype = anno$celltype,
                             border = TRUE, 
                             gap = unit(1,"mm"),
                             col = list(celltype = metadata(spe)$color_vectors$celltype),
                             which = "row")

ha_meta <- HeatmapAnnotation(n_cells = anno_barplot(anno$ncells, width = unit(10, "mm")),
                             n_PID = anno_barplot(n_PID, width = unit(10, "mm")),
                             indication = anno_barplot(indication,width = unit(10, "mm"),
                                                       gp = gpar(fill = metadata(spe)$color_vectors$indication)),
                             border = TRUE, 
                             annotation_name_rot = 90,
                             gap = unit(1,"mm"),
                             which = "row")

### 2.2 Spatial features

# Add number of neighbors to spe object (saved in colPair)
n_neighbors <- colPair(spe) %>% 
  as.data.frame() %>% 
  group_by(from) %>% 
  count() %>% 
  arrange(desc(n))

spe$n_neighbors <- n_neighbors$n[match(seq_along(colnames(spe)), n_neighbors$from)]
spe$n_neighbors <- spe$n_neighbors %>% replace_na(0)

# Select spatial features and average over celltypes
spatial <- colData(spe) %>% 
  as.data.frame() %>% 
  select(area, celltype, n_neighbors)

spatial <- spatial %>% 
  select(-celltype) %>% 
  aggregate(by = list(celltype = spatial$celltype), FUN = mean) %>% 
  column_to_rownames("celltype")

# Create HeatmapAnnotation object
ha_spatial <- HeatmapAnnotation(
  area = spatial$area,
  n_neighbors = spatial$n_neighbors,
  border = TRUE,
  gap = unit(1,"mm"),
  which = "row")

### 3. Plot rich heatmap ###

# Create HeatmapList object
h_list <- h_type +
  h_state +
  ha_anno +
  ha_spatial +
  ha_meta

# Add customized legend for anno_barplot()
lgd <- Legend(title = "indication", at = colnames(indication), 
              legend_gp = gpar(fill = metadata(spe)$color_vectors$indication))

# Plot
draw(h_list,annotation_legend_list = list(lgd))

#Finally, we save the updated SpatialExperiment object.
saveRDS(spe, "spe.rds")
saveRDS(spe_cat, "spe_cat.rds")


































