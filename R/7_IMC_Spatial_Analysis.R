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
library(SpatialExperiment)
library(tidyverse)
library(ggridges)

#1- Load data

#Highly multiplexed imaging technologies acquire the spatial distributions of molecule abundances across tissue sections. 
#As such, analyzing single cells in their spatial tissue context is a key strength of these technologies.
#A number of software packages such as squidpy, giotto and Seurat have been developed to analyse and visualize cells in their spatial context. 
#The following chapter will highlight the use of imcRtools and other Bioconductor tools tools to visualize and analyse single-cell data obtained from highly multiplexed imaging technologies.

#First, we will read in the previously generated SpatialExperiment object.
spe <- readRDS("spe.rds")

#2- Spatial interaction graphs

#The imcRtools package further allows the ad hoc consctruction of spatial graphs directly using a SpatialExperiment or SingleCellExperiment object 
#while considering the spatial location (centroids) of individual cells. 
#The buildSpatialGraph function allows constructing spatial graphs by detecting the k-nearest neighbors in 2D (knn), 
#by detecting all cells within a given distance to the center cell (expansion - similar to the radius option in Squidpy) and by Delaunay triangulation (delaunay).

#When constructing a knn graph, the number of neighbors (k) needs to be set and (optionally) the maximum distance to consider (max_dist) can be specified. 
#When constructing a graph via expansion, the distance to expand (threshold) needs to be provided. 
#For graphs constructed via Delaunay triangulation, the max_dist parameter can be set to avoid unusually large connections at the edge of the image.

library(imcRtools)
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "knn", k = 20)
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "expansion", threshold = 20)
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "delaunay", max_dist = 50)

#The spatial graphs are stored in colPair(spe, name) slots. 
#These slots store SelfHits objects representing edge lists 
#in which the first column indicates the index of the “from” cell and the second column the index of the “to” cell. 
#Each edge list is newly constructed when subsetting the object.

#Here, colPair(spe, "neighborhood") stores the spatial graph constructed by steinbock, 
#colPair(spe, "knn_interaction_graph") stores the knn spatial graph, 
#colPair(spe, "expansion_interaction_graph") stores the expansion graph and 
#colPair(spe, "delaunay_interaction_graph") stores the graph constructed by Delaunay triangulation.
colPairNames(spe)

#3- Spatial visualization

#Section 11 highlights the use of the cytomapper package to visualize multichannel images and segmentation masks. 
#Here, we introduce the plotSpatial function of the imcRtools package to visualize the cells’ centroids and cell-cell interactions as spatial graphs.

#In the following example, we select one image for visualization purposes. 
#Here, each dot (node) represents a cell and edges are drawn between cells in close physical proximity as detected by steinbock or the buildSpatialGraph function. 
#Nodes are variably colored based on the cell type and edges are colored in grey.
library(ggplot2)
library(viridis)

# steinbock interaction graph 
plotSpatial(spe[,spe$sample_id == "Patient3_001"], 
            node_color_by = "celltype", 
            img_id = "sample_id", 
            draw_edges = TRUE, 
            colPairName = "neighborhood", 
            nodes_first = FALSE, 
            edge_color_fix = "grey") + 
  scale_color_manual(values = metadata(spe)$color_vectors$celltype) +
  ggtitle("steinbock interaction graph")

# knn interaction graph 
plotSpatial(spe[,spe$sample_id == "Patient3_001"], 
            node_color_by = "celltype", 
            img_id = "sample_id", 
            draw_edges = TRUE, 
            colPairName = "knn_interaction_graph", 
            nodes_first = FALSE,
            edge_color_fix = "grey") + 
  scale_color_manual(values = metadata(spe)$color_vectors$celltype) +
  ggtitle("knn interaction graph")

# expansion interaction graph 
plotSpatial(spe[,spe$sample_id == "Patient3_001"], 
            node_color_by = "celltype", 
            img_id = "sample_id", 
            draw_edges = TRUE, 
            colPairName = "expansion_interaction_graph", 
            nodes_first = FALSE, 
            directed = FALSE,
            edge_color_fix = "grey") + 
  scale_color_manual(values = metadata(spe)$color_vectors$celltype) +
  ggtitle("expansion interaction graph")

# delaunay interaction graph 
plotSpatial(spe[,spe$sample_id == "Patient3_001"], 
            node_color_by = "celltype", 
            img_id = "sample_id", 
            draw_edges = TRUE, 
            colPairName = "delaunay_interaction_graph", 
            nodes_first = FALSE,
            edge_color_fix = "grey") + 
  scale_color_manual(values = metadata(spe)$color_vectors$celltype) +
  ggtitle("delaunay interaction graph")

#Nodes can also be colored based on the cells’ expression levels (e.g., E-cadherin expression) 
#and their size can be adjusted (e.g., based on measured cell area).
plotSpatial(spe[,spe$sample_id == "Patient3_001"], 
            node_color_by = "Ecad", 
            assay_type = "exprs",
            img_id = "sample_id", 
            draw_edges = TRUE, 
            colPairName = "expansion_interaction_graph", 
            nodes_first = FALSE, 
            node_size_by = "area", 
            directed = FALSE,
            edge_color_fix = "grey") + 
  scale_size_continuous(range = c(0.1, 2)) +
  ggtitle("E-cadherin expression")

#Finally, the plotSpatial function allows displaying all images at once. 
#This visualization can be useful to quickly detect larger structures of interest.
plotSpatial(spe, 
            node_color_by = "celltype", 
            img_id = "sample_id", 
            node_size_fix = 0.5) + 
  scale_color_manual(values = metadata(spe)$color_vectors$celltype)

#4- Spatial Community Analysis

#The detection of spatial communities was proposed by (Jackson et al. 2020). 
#Here, cells are clustered solely based on their interactions as defined by the spatial object graph. 
#In the following example, we perform spatial community detection separately for tumor and stromal cells.

#The general procedure is as follows:
#1. create a colData(spe) entry that specifies if a cell is part of the tumor or stroma compartment. 
#2. use the detectCommunity function of the imcRtools package to cluster cells within the tumor or stroma compartment solely based on their spatial interaction graph as constructed by the steinbock package.
#Both tumor and stromal spatial communities are stored in the colData of the SpatialExperiment object under the spatial_community identifier.
#We set the seed argument within the SerialParam function for reproducibility purposes. 
#This is important as the global seed is not recognized by functions provided by the BiocParallel package.

spe$tumor_stroma <- ifelse(spe$celltype == "Tumor", "Tumor", "Stroma")

library(BiocParallel)
spe <- detectCommunity(spe, 
                       colPairName = "neighborhood", 
                       size_threshold = 10,
                       group_by = "tumor_stroma",
                       BPPARAM = SerialParam(RNGseed = 220819))
#We can now separately visualize the tumor and stromal communities.
#Spatial Tumor communities
plotSpatial(spe[,spe$celltype == "Tumor"], 
            node_color_by = "spatial_community", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
  theme(legend.position = "none") +
  ggtitle("Spatial tumor communities") +
  scale_color_manual(values = rev(colors()))

#Spatial Stromal communities
plotSpatial(spe[,spe$celltype != "Tumor"], 
            node_color_by = "spatial_community", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
  theme(legend.position = "none") +
  ggtitle("Spatial non-tumor communities") +
  scale_color_manual(values = rev(colors()))

#In the next step, the fraction of cell types within each spatial stromal community is displayed.
library(pheatmap)
library(viridis)

for_plot <- prop.table(table(spe[,spe$celltype != "Tumor"]$spatial_community, spe[,spe$celltype != "Tumor"]$celltype), margin = 1)

pheatmap(for_plot, color = viridis(100), show_rownames = FALSE)

#5- Cellular neighborhood analysis

#The following section highlights the use of the imcRtools package to detect cellular neighborhoods. 
#This approach has been proposed by (Goltsev et al. 2018) and (Schürch et al. 2020) to group cells based on information contained in their direct neighborhood.
#(Goltsev et al. 2018) perfomed Delaunay triangulation-based graph construction, neighborhood aggregation and then clustered cells.
#(Schürch et al. 2020) on the other hand constructed a 10-nearest neighbor graph before aggregating information across neighboring cells.
#In the following code chunk we will use the 20-nearest neighbor graph as constructed above to define the direct cellular neighborhood. 
#The aggregateNeighbors function allows neighborhood aggregation in 2 different ways:

#A- For each cell the function computes the fraction of cells of a certain type (e.g., cell type) among its neighbors.
#B- For each cell it aggregates (e.g., mean) the expression counts across all neighboring cells.

#Based on these measures, cells can now be clustered into cellular neighborhoods. 
#We will first compute the fraction of the different cell types among the 20-nearest neighbors and use kmeans clustering to group cells into 6 cellular neighborhoods.

# By celltypes
spe <- aggregateNeighbors(spe, colPairName = "knn_interaction_graph", 
                          aggregate_by = "metadata", count_by = "celltype")

set.seed(220705)

cn_1 <- kmeans(spe$aggregatedNeighbors, centers = 6)
spe$cn_celltypes <- as.factor(cn_1$cluster)

plotSpatial(spe, 
            node_color_by = "cn_celltypes", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
  scale_color_brewer(palette = "Set3")

#The next code chunk visualizes the cell type compositions of the detected cellular neighborhoods (CN).
library(tidyverse)
for_plot <- colData(spe) %>% as_tibble() %>%
  group_by(cn_celltypes, celltype) %>%
  summarize(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  pivot_wider(id_cols = cn_celltypes, names_from = celltype, 
              values_from = freq, values_fill = 0) %>%
  ungroup() %>%
  select(-cn_celltypes)

pheatmap(for_plot, color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column")

#We will now detect cellular neighborhoods by computing the mean expression across the 20-nearest neighbor prior to kmeans clustering (k=6).
# By expression
spe <- aggregateNeighbors(spe, colPairName = "knn_interaction_graph", 
                          aggregate_by = "expression", assay_type = "exprs",
                          subset_row = rowData(spe)$use_channel)
cn_2 <- kmeans(spe$mean_aggregatedExpression, centers = 6)
spe$cn_expression <- as.factor(cn_2$cluster)

plotSpatial(spe, 
            node_color_by = "cn_expression", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
  scale_color_brewer(palette = "Set3")

#we can visualize the cell type composition of each cellular neighborhood.
for_plot <- colData(spe) %>% as_tibble() %>%
  group_by(cn_expression, celltype) %>%
  summarize(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  pivot_wider(id_cols = cn_expression, names_from = celltype, 
              values_from = freq, values_fill = 0) %>%
  ungroup() %>%
  select(-cn_expression)

pheatmap(for_plot, color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column")

#When clustering cells based on the mean expression within the direct neighborhood, 
#tumor patches are split across CN 4 and CN 6 without forming a clear tumor/stroma interface. 
#This result reflects patient-to-patient differences in the expression of tumor markers.
#Of note: constructing a 20-nearest neighbor graph and clustering using kmeans with k=6 is only an example. 
#Similar to the analysis done in Section 9.2.2, it is recommended to perform a parameter sweep across 
#different graph construction algorithms and different parmaters k for kmeans clustering. 
#Finding the best CN detection settings is also subject to the question at hand. 
#Constructing graphs with more neighbors usually results in larger CNs.

#6- lisaClust
#An alternative to the aggregateNeighbors function is provided by the lisaClust Bioconductor package (Patrick et al. 2021). 
#In contrast to imcRtools, the lisaClust package computes local indicators of spatial associations (LISA) functions and clusters cells based on those. 
#More precise, the package summarizes L-functions from a Poisson point process model to derive numeric vectors 
#for each cell which can then again be clustered using kmeans.

#The lisa function requires a SegmentedCells object which can be generated using the spicyR package.
library(lisaClust)
library(spicyR)

cells <- data.frame(row.names = colnames(spe))
cells$ObjectNumber <- spe$ObjectNumber
cells$ImageNumber <- spe$sample_id
cells$AreaShape_Center_X <- spatialCoords(spe)[,"Pos_X"]
cells$AreaShape_Center_Y <- spatialCoords(spe)[,"Pos_Y"]
cells$cellType <- spe$celltype

lisa_sc <- SegmentedCells(cells, cellProfiler = TRUE)
lisa_sc

#After creating the SegmentedCells object, the lisa function computes LISA curves across a given set of distances. 
#In the following example, we calculate the LISA curves within a 10µm, 20µm and 50µm neighborhood around each cell. 
#Increasing these radii will lead to broader and smoother spatial clusters.
#However, a number of parameter settings should be tested to estimate the robustness of the results.
lisaCurves <- lisa(lisa_sc, Rs = c(10, 20, 50))

# Set NA to 0
lisaCurves[is.na(lisaCurves)] <- 0

lisa_clusters <- kmeans(lisaCurves, centers = 6)$cluster

spe$lisa_clusters <- as.factor(lisa_clusters)

plotSpatial(spe, 
            node_color_by = "lisa_clusters", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
  scale_color_brewer(palette = "Set3")

#We can now observe the cell type composition per spatial cluster.
for_plot <- colData(spe) %>% as_tibble() %>%
  group_by(lisa_clusters, celltype) %>%
  summarize(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  pivot_wider(id_cols = lisa_clusters, names_from = celltype, 
              values_from = freq, values_fill = 0) %>%
  ungroup() %>%
  select(-lisa_clusters)

pheatmap(for_plot, color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column")

#7- Spatial Context Analysis

#Downstream of CN assignments, we will analyze the spatial context (SC) of each cell using three functions from imcRtools.
#While CNs can represent sites of unique local processes, the term SC was coined by Bhate and colleagues (Bhate et al. 2022) 
#and describes tissue regions in which distinct CNs may be interacting. 
#Hence, SCs may be interesting regions of specialized biological events.

#Here, we will first detect SCs using the detectSpatialContext function. 
#This function relies on CN fractions for each cell in a spatial interaction graph (originally a KNN graph),
#which we will calculate using buildSpatialGraph and aggregateNeighbors. 
#We will focus on the CNs derived from cell type fractions but other CN assignments are possible.

#Of note, the window size (k for KNN) for buildSpatialGraph should reflect a length scale 
#on which biological signals can be exchanged and depends, among others, on cell density and tissue area. 
#In view of their divergent functionality, we recommend to use a larger window size for SC (interaction between local processes) than for CN (local processes) detection. 
#Since we used a 20-nearest neighbor graph for CN assignment, we will use a 40-nearest neighbor graph for SC detection. 
#As before, different parameters should be tested.
#Subsequently, the CN fractions are sorted from high-to-low and the SC of each cell is assigned as the minimal combination of SCs that additively surpass a user-defined threshold. 
#The default threshold of 0.9 aims to represent the dominant CNs, hence the most prevalent signals, in a given window.
#For more details and biological validation, please refer to (Bhate et al. 2022).

library(circlize)
library(RColorBrewer)

# Generate k-nearest neighbor graph for SC detection (k=40) 
spe <- buildSpatialGraph(spe, img_id = "sample_id", 
                         type = "knn", 
                         name = "knn_spatialcontext_graph", 
                         k = 40)

# Aggregate based on clustered_neighbors
spe <- aggregateNeighbors(spe, 
                          colPairName = "knn_spatialcontext_graph",
                          aggregate_by = "metadata",
                          count_by = "cn_celltypes",
                          name = "aggregatedNeighborhood")

# Detect spatial contexts
spe <- detectSpatialContext(spe, 
                            entry = "aggregatedNeighborhood",
                            threshold = 0.90,
                            name = "spatial_context")

# Define SC color scheme
col_SC <- setNames(colorRampPalette(brewer.pal(9, "Paired"))(length(unique(spe$spatial_context))), 
                   sort(unique(spe$spatial_context)))

# Visualize spatial contexts on images
plotSpatial(spe, 
            node_color_by = "spatial_context", 
            img_id = "sample_id", 
            node_size_fix = 0.5, 
            colPairName = "knn_spatialcontext_graph") +
  scale_color_manual(values = col_SC)

#For ease of interpretation, we will directly compare the CN and SC assignments for Patient3_001.
library(patchwork)

# Compare CN and SC for one patient 
#As expected, we can observe that interfaces between different CNs make up distinct SCs. 
#For instance, interface between CN 3 (TLS region consisting of B and BnT cells) and CN 4 (Plasma- and T-cell dominated) turns to SC 3_4. 
#On the other hand, the core of CN 3 becomes SC 3, since for the neighborhood for these cells is just the cellular neighborhood itself.
p1 <- plotSpatial(spe[,spe$sample_id == "Patient3_001"], 
                  node_color_by = "cn_celltypes", 
                  img_id = "sample_id", 
                  node_size_fix = 0.5, 
                  colPairName = "knn_interaction_graph") +
  scale_color_brewer(palette = "Set3")

p2 <- plotSpatial(spe[,spe$sample_id == "Patient3_001"], 
                  node_color_by = "spatial_context", 
                  img_id = "sample_id", 
                  node_size_fix = 0.5, 
                  colPairName = "knn_spatialcontext_graph") +
  scale_color_manual(values = col_SC, limits = force)

p1 + p2

#Next, we filter the SCs based on user-defined thresholds for number of group entries (here at least 3 patients) 
#and/or total number of cells (here minimum of 100 cells) per SC with filterSpatialContext.

## Filter spatial contexts
# By number of group entries
spe <- filterSpatialContext(spe, 
                            entry = "spatial_context",
                            group_by = "patient_id", 
                            group_threshold = 3)

plotSpatial(spe, 
            node_color_by = "spatial_context_filtered", 
            img_id = "sample_id", 
            node_size_fix = 0.5, 
            colPairName = "knn_spatialcontext_graph") +
  scale_color_manual(values = col_SC, limits = force)

# By number of group entries and total number of cells
spe <- filterSpatialContext(spe, 
                            entry = "spatial_context",
                            group_by = "patient_id", 
                            group_threshold = 3,
                            cells_threshold = 100)

plotSpatial(spe, 
            node_color_by = "spatial_context_filtered", 
            img_id = "sample_id", 
            node_size_fix = 0.5, 
            colPairName = "knn_spatialcontext_graph") +
  scale_color_manual(values = col_SC, limits = force)

#Lastly, we can use the plotSpatialContext function to generate SC graphs, 
#analogous to CN combination maps in (Bhate et al. 2022). 
#Returned objects are ggplots, which can be easily modified further. 
#We will create a SC graph for the filtered SCs here.

## Plot spatial context graph 

# Colored by name and size by n_cells
plotSpatialContext(spe, 
                   entry = "spatial_context_filtered",
                   group_by = "sample_id",
                   node_color_by = "name",
                   node_size_by = "n_cells",
                   node_label_color_by = "name")


# Colored by n_cells and size by n_group                   
plotSpatialContext(spe, 
                   entry = "spatial_context_filtered",
                   group_by = "sample_id",
                   node_color_by = "n_cells",
                   node_size_by = "n_group",
                   node_label_color_by = "n_cells") +
  scale_color_viridis()


#8-  Patch Detection Analysis

#The previous section focused on detecting cellular neighborhoods in a rather unsupervised fashion. 
#However, the imcRtools package also provides methods for detecting spatial compartments in a supervised fashion. 
#The patchDetection function allows the detection of connected sets of similar cells as proposed by (Hoch et al. 2022). 
#In the following example, we will use the patchDetection function to detect function to detect tumor patches in three steps:

#A- Find connected sets of tumor cells (using the steinbock graph).
#B- Components which contain less than 10 cells are excluded.
#C- Expand the components by 1µm to construct a concave hull around the patch and include cells within the patch.

spe <- patchDetection(spe, 
                      patch_cells = spe$celltype == "Tumor",
                      img_id = "sample_id",
                      expand_by = 1,
                      min_patch_size = 10,
                      colPairName = "neighborhood")

plotSpatial(spe, 
            node_color_by = "patch_id", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
  theme(legend.position = "none") +
  scale_color_manual(values = colors())

#We can now calculate the fraction of T cells within each tumor patch to roughly estimate T cell infiltration.
library(tidyverse)
colData(spe) %>% as_tibble() %>%
  group_by(patch_id, sample_id) %>%
  summarize(Tcell_count = sum(celltype == "CD8 T cell" | celltype == "CD4 T cell"),
            patch_size = n(),
            Tcell_freq = Tcell_count / patch_size) %>%
  ggplot() +
  geom_point(aes(log10(patch_size), Tcell_freq, color = sample_id)) +
  theme_classic()

#We can now measure the size of each patch using the patchSize function and visualize tumor patch distribution per patient.
patch_size <- patchSize(spe, "patch_id")

patch_size <- merge(patch_size, 
                    colData(spe)[match(patch_size$patch_id, spe$patch_id),], 
                    by = "patch_id")

ggplot(as.data.frame(patch_size)) + 
  geom_boxplot(aes(patient_id, log10(size))) +
  geom_point(aes(patient_id, log10(size)))


#The minDistToCells function can be used to calculate the minimum distance between each cell and a cell set of interest. 
#Here, we highlight its use to calculate the minimum distance of all cells to the detected tumor patches. 
#Negative values indicate the minimum distance of each tumor patch cell to a non-tumor patch cell.
spe <- minDistToCells(spe, 
                      x_cells = !is.na(spe$patch_id), 
                      img_id = "sample_id")

plotSpatial(spe, 
            node_color_by = "distToCells", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
  scale_color_gradient2(low = "dark blue", mid = "white", high = "dark red")

#Finally, we can observe the minimum distances to tumor patches in a cell type specific manner.
library(ggridges)
ggplot(as.data.frame(colData(spe))) + 
  geom_density_ridges(aes(distToCells, celltype, fill = celltype)) +
  geom_vline(xintercept = 0, color = "dark red", size = 2) +
  scale_fill_manual(values = metadata(spe)$color_vectors$celltype)

#9- Interaction analysis

#The next section focuses on statistically testing the pairwise interaction between all cell types of the dataset. 
#For this, the imcRtools package provides the testInteractions function which implements the interaction testing strategy proposed by (Schapiro et al. 2017).
#Per grouping level (e.g., image), the testInteractions function computes the averaged cell type/cell type interaction count 
#and computes this count against an empirical null distribution which is generated by permuting all cell labels (while maintaining the tissue structure).
#In the following example, we use the steinbock generated spatial interaction graph and estimate the interaction or avoidance between cell types in the dataset.

library(scales)
out <- testInteractions(spe, 
                        group_by = "sample_id",
                        label = "celltype", 
                        colPairName = "neighborhood",
                        BPPARAM = SerialParam(RNGseed = 221029))

head(out)

#The returned DataFrame contains the test results per grouping level (in this case the image ID, group_by),
#“from” cell type (from_label) and “to” cell type (to_label). 
#The sigval entry indicates if a pair of cell types is significantly interacting (sigval = 1),
#if a pair of cell types is significantly avoiding (sigval = -1) or if no significant interaction or avoidance was detected.
#These results can be visualized by computing the sum of the sigval entries across all images:
out %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#The imcRtools package further implements an interaction testing strategy proposed by (Schulz et al. 2018) 
#where the hypothesis is tested if at least n cells of a certain type are located around a target cell type (from_cell). 
#This type of testing can be performed by selecting method = "patch" and specifying the number of patch cells via the patch_size parameter.

out <- testInteractions(spe, 
                        group_by = "sample_id",
                        label = "celltype", 
                        colPairName = "neighborhood",
                        method = "patch", 
                        patch_size = 3,
                        BPPARAM = SerialParam(RNGseed = 221029))

out %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#These results are comparable to the interaction testing presented above. 
#The main difference comes from the lack of symmetry. We can now for example see that 3 or more myeloid cells sit around CD4 T cells 
#while this interaction is not as strong when considering CD4 T cells sitting around myeloid cells.

#Finally, we save the updated SpatialExperiment object.
saveRDS(spe, "spe.rds")


















































































































