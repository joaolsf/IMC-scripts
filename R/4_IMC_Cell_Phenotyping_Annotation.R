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
library(garnett)
library(monocle3)
library(Seurat)
library(SeuratDisk)
library(basilisk)
library(zellkonverter)
library(reticulate)
library(monocle)
library(Matrix)
library(org.Hs.eg.db)


#1- Load data
#First, we will read in the previously generated SpatialExperiment object.
spe <- readRDS("spe.rds")
spe

#2- Cell phenotyping
#2.1- Sample cells
set.seed(220619)
cur_cells <- sample(seq_len(ncol(spe)), 2000)

# 2.2- Clustering approaches
# 2.2.1- Rphenograph
#we select the asinh-transformed mean pixel intensities per cell and channel and 
#subset the channels to the ones containing biological variation. 
#This matrix is transposed to store cells in rows. 
#Within the Rphenograph function, we select the 45 nearest neighbors for graph building and 
#louvain community detection (default). 
#The function returns a list of length 2, the first entry being the graph and the second entry containing the community object. 
#Calling membership on the community object will return cluster IDs for each cell. 
#These cluster IDs are then stored within the colData of the SpatialExperiment object. 
#Cluster IDs are mapped on top of the UMAP embedding and single-cell marker expression within each cluster are visualized in form of a heatmap.
library(Rphenograph)
library(igraph)
library(dittoSeq)
library(viridis)

mat <- t(assay(spe, "exprs")[rowData(spe)$use_channel,])

out <- Rphenograph(mat, k = 45)

clusters <- factor(membership(out[[2]]))

spe$pg_clusters <- clusters

dittoDimPlot(spe, var = "pg_clusters", 
             reduction.use = "UMAP", size = 0.2,
             do.label = TRUE) +
  ggtitle("Phenograph clusters expression on UMAP")

#heatmap of the Phenograph clusters
dittoHeatmap(spe[,cur_cells], 
             genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "exprs", scale = "none",
             heatmap.colors = viridis(100), 
             annot.by = c("pg_clusters", "patient_id"),
             annot.colors = c(dittoColors(1)[1:length(unique(spe$pg_clusters))],
                              metadata(spe)$color_vectors$patient_id))

#We can observe that some of the clusters only contain cells of a single patient. 
#This can often be observed in the tumor compartment.
#In the next step, we use the integrated cells (see Section 8) in low dimensional embedding for clustering. 
#Here, the low dimensional embedding can be directly accessed from the reducedDim slot.
mat <- reducedDim(spe, "fastMNN")

out <- Rphenograph(mat, k = 45)

clusters <- factor(membership(out[[2]]))

spe$pg_clusters_corrected <- clusters

dittoDimPlot(spe, var = "pg_clusters_corrected", 
             reduction.use = "UMAP_mnnCorrected", size = 0.2,
             do.label = TRUE) +
  ggtitle("Phenograph clusters expression on UMAP, integrated cells")

#heatmao of the Phenograph clusters
dittoHeatmap(spe[,cur_cells], 
             genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "exprs", scale = "none",
             heatmap.colors = viridis(100),
             annot.by = c("pg_clusters_corrected","patient_id"), #"pg_clusters",
             annot.colors = c(dittoColors(1)[1:length(unique(spe$pg_clusters_corrected))],
                              #dittoColors(1)[1:length(unique(spe$pg_clusters))],
                              metadata(spe)$color_vectors$patient_id))

# 2.2.2- Shared nearest neighbour graph
library(bluster)
library(BiocParallel)
library(ggplot2)
#we sample 10000 cells but it is recommended to increase this number to include rare cells into parameter testing.
sam_cells <- sample(seq_len(ncol(spe)), 10000)

#We test two different settings for k, two for type and two for cluster.fun. 
#This function call is parallelized by setting the BPPARAM parameter. 
#We use the approxSilhouette function to compute the silhouette width for each cell and compute the average across all cells per parameter setting.
#Please see ?silhouette for more information on how the silhouette width id computed for each cell. 
#A large average silhouette width indicates cells that are well clustered.
mat <- t(assay(spe, "exprs")[rowData(spe)$use_channel,])

combinations <- clusterSweep(mat[sam_cells,], BLUSPARAM=SNNGraphParam(),
                             k=c(10L, 20L), 
                             type = c("rank", "jaccard"), 
                             cluster.fun=c("walktrap", "louvain"),
                             BPPARAM = MulticoreParam(RNGseed = 220427))

sil <- vapply(as.list(combinations$clusters), 
              function(x) mean(approxSilhouette(mat[sam_cells,], x)$width), 
              0)

ggplot(data.frame(method = names(sil),
                  sil = sil)) +
  geom_point(aes(method, sil)) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Cluster parameter combination") +
  ylab("Average silhouette width")

#Once parameter settings are known, we can either use the clusterRows function of the bluster package to cluster cells or 
#its convenient wrapper function exported by the scran package. 
#The scran::clusterCells function accepts a SpatialExperiment (or SingleCellExperiment) object which stores cells in columns. 
#By default, the function detects the 10 nearest neighbours for each cell, performs rank-based weighting of edges (see ?makeSNNGraph for more information) 
#and uses the cluster_walktrap function to detect communities in the graph.

#Running clustering in the non-batch corrected reduction
library(scran)

clusters <- clusterCells(spe[rowData(spe)$use_channel,], 
                         assay.type = "exprs", 
                         BLUSPARAM = SNNGraphParam(k=20, 
                                                   cluster.fun = "louvain",
                                                   type = "rank",
                                                   BPPARAM = bpparam()))

spe$nn_clusters <- clusters

dittoDimPlot(spe, var = "nn_clusters", 
             reduction.use = "UMAP", size = 0.2,
             do.label = TRUE) +
  ggtitle("SNN clusters expression on UMAP")

dittoHeatmap(spe[,cur_cells], 
             genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "exprs", scale = "none",
             heatmap.colors = viridis(100), 
             annot.by = c("nn_clusters", "patient_id"),
             annot.colors = c(dittoColors(1)[1:length(unique(spe$nn_clusters))],
                              metadata(spe)$color_vectors$patient_id))


#This function was used by (Tietscher et al. 2022) to cluster cells obtained by IMC. 
#Setting type = "jaccard" performs clustering similar to Rphenograph above and Seurat.
#Similar to the results obtained by Rphenograph, some of the clusters are patient-specific. 
#We can now perform clustering of the integrated cells by directly specifying which low-dimensional embedding to use:
clusters <- clusterCells(spe[rowData(spe)$use_channel,], 
                         use.dimred = "fastMNN", 
                         BLUSPARAM = NNGraphParam(k=20, 
                                                  cluster.fun = "louvain",
                                                  type = "rank",
                                                  BPPARAM = bpparam()))

spe$nn_clusters_corrected <- clusters

#UMAP embedding of the clusters of the integrated (batch-corrected) cells
dittoDimPlot(spe, var = "nn_clusters_corrected", 
             reduction.use = "UMAP_mnnCorrected", size = 0.2,
             do.label = TRUE) +
  ggtitle("SNN clusters expression on UMAP, integrated cells")

#Heatmap of the clusters of the integrated (batch-corrected) cells
dittoHeatmap(spe[,cur_cells], 
             genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "exprs", scaled.to.max = 'FALSE', scale = 'none',
             heatmap.colors = viridis(100), 
             annot.by = c("nn_clusters_corrected", "nn_clusters","patient_id"),
             annot.colors = c(dittoColors(1)[1:length(unique(spe$nn_clusters_corrected))],
                              dittoColors(1)[1:length(unique(spe$nn_clusters))],
                              metadata(spe)$color_vectors$patient_id))

# 2.2.3- Self organizing maps
#An alternative to graph-based clustering is offered by the CATALYST package. 
#The cluster function internally uses the FlowSOM package to group cells into 100 (default) clusters based on self organizing maps (SOM).
#In the next step, the ConsensusClusterPlus package is used to perform hierarchical consensus clustering of the previously detected 100 SOM nodes into 1 to maxK clusters. 
#Cluster stability for each k can be assessed by plotting the delta_area(spe). 
#The optimal number of clusters can be found by selecting the k at which a plateau is reached. 
#In the example below, the optimal k lies somewhere around 13.
library(CATALYST)

# Run FlowSOM and ConsensusClusterPlus clustering
spe <- cluster(spe, 
               features = rownames(spe)[rowData(spe)$use_channel],
               maxK = 30,
               seed = 220410)

# Assess cluster stability
delta_area(spe)

#Running clustering in the non-batch corrected reduction
spe$som_clusters <- cluster_ids(spe, "meta13")

dittoDimPlot(spe, var = "som_clusters", 
             reduction.use = "UMAP", size = 0.2,
             do.label = TRUE) +
  ggtitle("SOM clusters expression on UMAP")

dittoHeatmap(spe[,cur_cells], 
             genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "exprs", scale = "none",
             heatmap.colors = viridis(100), 
             annot.by = c("som_clusters", "patient_id"),
             annot.colors = c(dittoColors(1)[1:length(unique(spe$som_clusters))],
                              metadata(spe)$color_vectors$patient_id))


#The CATALYST package does not provide functionality to perform FlowSOM and ConsensusClusterPlus clustering directly on the batch-corrected, integrated cells. 
#As an alternative to the CATALYST package, the bluster package provides SOM clustering when specifying the SomParam() parameter. 
#Similar to the CATALYST approach, we will first cluster the dataset into 100 clusters (also called “codes”). 
#These codes are then further clustered into a maximum of 30 clusters using ConsensusClusterPlus (using hierarchical clustering and euclidean distance). 
#The delta area plot can be accessed using the (not exported) .plot_delta_area function from CATALYST.
#Here, it seems that the plateau is reached at a k of 16 and we will store the final cluster IDs within the SpatialExperiment object.
library(kohonen)
library(ConsensusClusterPlus)

# Select integrated cells
mat <- reducedDim(spe, "fastMNN")

# Perform SOM clustering
som.out <- clusterRows(mat, SomParam(100), full = TRUE)

# Cluster the 100 SOM codes into larger clusters
ccp <- ConsensusClusterPlus(t(som.out$objects$som$codes[[1]]),
                            maxK = 30,
                            reps = 100, 
                            distance = "euclidean", 
                            seed = 220410, 
                            plot = NULL)

# Visualize delta area plot
CATALYST:::.plot_delta_area(ccp)

# Link ConsensusClusterPlus clusters with SOM codes and save in object
#The FlowSOM clustering approach has been used by (Hoch et al. 2022) to sub-cluster tumor cells as measured by IMC.
som.cluster <- ccp[[16]][["consensusClass"]][som.out$clusters]
spe$som_clusters_corrected <- as.factor(som.cluster)

dittoDimPlot(spe, var = "som_clusters_corrected", 
             reduction.use = "UMAP_mnnCorrected", size = 0.2,
             do.label = TRUE) +
  ggtitle("Phenograph clusters expression on UMAP, integrated cells")

dittoHeatmap(spe[,cur_cells], 
             genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "exprs", scale = "none",
             heatmap.colors = viridis(100), 
             annot.by = c("som_clusters_corrected", "som_clusters","patient_id"),
             annot.colors = c(dittoColors(1)[1:length(unique(spe$som_clusters_corrected))],
                              dittoColors(1)[1:length(unique(spe$som_clusters))],
                              metadata(spe)$color_vectors$patient_id))

# 2.2.4- Compare between clustering approaches

#Finally, we can compare the results of different clustering approaches.
#For this, we visualize the number of cells that are shared between different clustering results in a pairwise fashion. 
#In the following heatmaps a high match between clustering results can be seen for those clusters 
#that are uniquely detected in both approaches.

#First, we will visualize the match between the three different approaches applied to the asinh-transformed counts.
library(patchwork)
library(pheatmap)
library(gridExtra)

tab1 <- table(paste("Rphenograph", spe$pg_clusters), 
              paste("SNN", spe$nn_clusters))
tab2 <- table(paste("Rphenograph", spe$pg_clusters), 
              paste("SOM", spe$som_clusters))
tab3 <- table(paste("SNN", spe$nn_clusters), 
              paste("SOM", spe$som_clusters))

pheatmap(log10(tab1 + 10), color = viridis(100))
pheatmap(log10(tab2 + 10), color = viridis(100))
pheatmap(log10(tab3 + 10), color = viridis(100))

#Below, a comparison between the clustering results of the integrated cells is shown.
#In comparison to clustering on the non-integrated cells, the clustering results of the integrated cells show higher overlap. 
#The SNN approach resulted in fewer clusters and therefore matches better with the SOM clustering approach.
tab1 <- table(paste("Rphenograph", spe$pg_clusters_corrected), 
              paste("SNN", spe$nn_clusters_corrected))
tab2 <- table(paste("Rphenograph", spe$pg_clusters_corrected), 
              paste("SOM", spe$som_clusters_corrected))
tab3 <- table(paste("SNN", spe$nn_clusters_corrected), 
              paste("SOM", spe$som_clusters_corrected))

pheatmap(log10(tab1 + 10), color = viridis(100))
pheatmap(log10(tab2 + 10), color = viridis(100))
pheatmap(log10(tab3 + 10), color = viridis(100))


#Finally, we save the updated SpatialExperiment object.
saveRDS(spe, "spe.rds")



#3- Cell type classification

#In this section, we will highlight a cell type classification approach based on ground truth labeling and random forest classification. 
#The rational for this supervised cell phenotyping approach is to use the information contained in the pre-defined markers to detect cells of interest.
#As related approaches, Astir and Garnett use pre-defined panel information to classify cell phenotypes based on their marker expression.

#3.1- Manual labeling of cells

#The cytomapper package provides the cytomapperShiny function that allows gating of cells based on their marker expression 
#and visualization of selected cells directly on the images.
library(cytomapper)
if (interactive()) {
  
  images <- readRDS("images.rds")
  masks <- readRDS("masks.rds")
  
  cytomapperShiny(object = spe, mask = masks, image = images, 
                  cell_id = "ObjectNumber", img_id = "sample_id")
}

#3.1.1- Define color vectors
#For consistent visualization of cell types, we will now pre-define their colors:
celltype <- setNames(c("#3F1B03", "#F4AD31", "#894F36", "#1C750C", "#EF8ECC", 
                       "#6471E2", "#4DB23B", "grey", "#F4800C", "#BF0A3D", "#066970"),
                     c("Tumor", "Stroma", "Myeloid", "CD8", "Plasma_cell", 
                       "Treg", "CD4", "undefined", "BnTcell", "Bcell", "Neutrophil"))

metadata(spe)$color_vectors$celltype <- celltype

#3.1.2- Read in and consolidate labeled data
#Here, we will read in the individual SpatialExperiment objects containing the labeled cells and concatenate them. 
#In the process of concatenating the SpatialExperiment objects along their columns, 
#the sample_id entry is appended by .1, .2, .3, ... due to replicated entries.
library(SingleCellExperiment)
label_files <- list.files("gated_cells",
                          full.names = TRUE, pattern = ".rds$")
# Read in SPE objects
spes <- lapply(label_files, readRDS)

# Merge SPE objects
concat_spe <- do.call("cbind", spes)

#In the following code chunk we will identify cells that were labeled multiple times. 
#This occurs when different cell phenotypes are gated per image and can affect immune cells
#that are located inside the tumor compartment.

#We will first identify those cells that were uniquely labeled. 
#In the next step, we will identify those cells that were labeled twice AND were labeled as Tumor cells. 
#These cells will be assigned their immune cell label. 
#Finally, we will save the unique labels within the original SpatialExperiment object.

cur_tab <- unclass(table(colnames(concat_spe), 
                         concat_spe$cytomapper_CellLabel))
cur_labels <- rep("doublets", nrow(cur_tab))
names(cur_labels) <- rownames(cur_tab)

# Single assignments
single_index <- rowSums(cur_tab) == 1
cur_labels[single_index] <- colnames(cur_tab)[apply(cur_tab[single_index,], 1, 
                                                    which.max)]

# Double assignment within the tumor
double_index <- rowSums(cur_tab) == 2 & cur_tab[,"Tumor"] == 1
no_tumor <- cur_tab[,colnames(cur_tab) != "Tumor"]
cur_labels[double_index] <- colnames(no_tumor)[apply(no_tumor[double_index,], 1, 
                                                     which.max)]

# Remove doublets
cur_labels <- cur_labels[cur_labels != "doublets"]
table(cur_labels)

# Transfer labels to SPE object
spe_labels <- rep("unlabeled", ncol(spe))
names(spe_labels) <- colnames(spe)
spe_labels[names(cur_labels)] <- cur_labels
spe$cell_labels <- spe_labels

# Number of cells labeled per patient
table(spe$cell_labels, spe$patient_id)

#3.1.3- Train classifier
#In this section, we will use the caret framework for machine learning in R. 
#This package provides an interface to train a number of regression and classification models in a coherent fashion.
#We use a random forest classifier due to low number of parameters, high speed and an observed high performance for cell type classification (Hoch et al. 2022).

#In the following section, we will first split the SpatialExperiment object into labeled and unlabeled cells. 
#Based on the labeled cells, we split the data into a train (75% of the data) and test (25% of the data) dataset. 
#We currently do not provide an independently labeled validation dataset.

#The caret package provides the trainControl function, which specifies model training parameters and the train function, 
#which performs the actual model training. 
#While training the model, we also want to estimate the best model parameters. 
#In the case of the chosen random forest model (method = "rf"), we only need to estimate a single parameters (mtry) which corresponds to the number of variables randomly sampled as candidates at each split.
#To estimate the best parameter, we will perform a 5-fold cross validation (set within trainControl) over a tune length of 5 entries to mtry.
library(caret)

# Split between labeled and unlabeled cells
lab_spe <- spe[,spe$cell_labels != "unlabeled"]
unlab_spe <- spe[,spe$cell_labels == "unlabeled"]

# Randomly split into train and test data
set.seed(221029)
trainIndex <- createDataPartition(factor(lab_spe$cell_labels), p = 0.75)
train_spe <- lab_spe[,trainIndex$Resample1]
test_spe <- lab_spe[,-trainIndex$Resample1]

# Specify train parameters for 5-fold cross validation
fitControl <- trainControl(method = "cv",
                           number = 5)

# Select the data for training
cur_mat <- t(assay(train_spe, "exprs")[rowData(train_spe)$use_channel,])

# Train a random forest model for predicting cell labels
# This call also performs parameter tuning
rffit <- train(x = cur_mat, 
               y = factor(train_spe$cell_labels),
               method = "rf", ntree = 1000,
               tuneLength = 5,
               trControl = fitControl)

rffit

# Classifier performance
ggplot(rffit) + 
  geom_errorbar(data = rffit$results,
                aes(ymin = Accuracy - AccuracySD,
                    ymax = Accuracy + AccuracySD),
                width = 0.4) +
  theme_classic(base_size = 15)
plot(varImp(rffit))

cm <- confusionMatrix(data = cur_pred, 
                      reference = factor(test_spe$cell_labels), 
                      mode = "everything")

cm
#To easily visualize these results, we can now plot the true positive rate (sensitivity) versus the false positive rate (1 - specificity):
library(tidyverse)

data.frame(cm$byClass) %>%
  mutate(class = sub("Class: ", "", rownames(cm$byClass))) %>%
  ggplot() + 
  geom_point(aes(1 - Specificity, Sensitivity, 
                 size = Detection.Rate,
                 fill = class),
             shape = 21) + 
  scale_fill_manual(values = metadata(spe)$color_vectors$celltype) +
  theme_classic(base_size = 15) + 
  ylab("Sensitivity (TPR)") +
  xlab("1 - Specificity (FPR)")

#Finally, to observe which cell phenotypes were wrongly classified,
#we can visualize the distribution of classification probabilities per cell phenotype class:
cur_pred <- predict(rffit, 
                    newdata = cur_mat, 
                    type = "prob")
cur_pred$truth <- factor(test_spe$cell_labels)

cur_pred %>%
  pivot_longer(cols = Bcell:Tumor) %>%
  ggplot() +
  geom_boxplot(aes(x = name, y = value, fill = name), outlier.size = 0.5) +
  facet_wrap(. ~ truth, ncol = 1) + 
  scale_fill_manual(values = metadata(spe)$color_vectors$celltype)  +
  theme(panel.background = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))

#3.1.4- Classification of new data
#In the final section, we will now use the tuned and tested random forest classifier to predict the cell phenotypes of the unlabeled data.
#First, we predict the cell phenotypes and extract their classification probabilities.

# Select unlabeled data
cur_mat <- t(assay(unlab_spe, "exprs")[rowData(unlab_spe)$use_channel,])

# Predict cell phenotypes
cell_class <- as.character(predict.train(rffit, 
                                         newdata = cur_mat, 
                                         type = "raw"))
names(cell_class) <- rownames(cur_mat)

table(cell_class)

# Extract classification probabilities
cell_prob <- predict.train(rffit, 
                           newdata = cur_mat, 
                           type = "prob")

library(ggridges)

# Distribution of maximum probabilities
tibble(max_prob = rowMax(as.matrix(cell_prob)),
       type = cell_class) %>%
  ggplot() +
  geom_density_ridges(aes(x = max_prob, y = cell_class, fill = cell_class)) +
  scale_fill_manual(values = metadata(spe)$color_vectors$celltype) +
  theme_classic(base_size = 15) +
  xlab("Maximum probability") +
  ylab("Cell type") + 
  xlim(c(0,1.2))

# Label undefined cells
cell_class[rowMax(as.matrix(cell_prob)) < 0.4] <- "undefined"

# Store labels in SpatialExperiment onject
cell_labels <- spe$cell_labels
cell_labels[colnames(unlab_spe)] <- cell_class
spe$celltype <- cell_labels 

table(spe$celltype, spe$patient_id)

#We next compare the cell classification against clustering results using the integrated cells.
tab1 <- table(spe$celltype, 
              paste("Rphenograph", spe$pg_clusters_corrected))
tab2 <- table(spe$celltype, 
              paste("SNN", spe$nn_clusters_corrected))
tab3 <- table(spe$celltype, 
              paste("SOM", spe$som_clusters_corrected))

pheatmap(log10(tab1 + 10), color = viridis(100))



#3.2- Cell Annotation with Garnett

#3.2.1- Build the CDS object - see Monocle documentation for help
#http://cole-trapnell-lab.github.io/monocle-release/docs/#the-celldataset-class

#The CellDataSet class:
#Because Garnett builds on Monocle 3, data for Garnett is held in objects of the cell_data_set (CDS) class. 
#This class is derived from the Bioconductor SingleCellExperiment class, which provides a common interface familiar to those who have analyzed single-cell data with Bioconductor. 
#Monocle 3 provides detailed documentation about how to generate an input CDS here.

#Monocle holds single cell expression data in objects of the CellDataSet class. 
#The class requires three input files:
#exprs, a numeric matrix of expression values, where rows are genes, and columns are cells
#phenoData, an AnnotatedDataFrame object, where rows are cells, and columns are cell attributes (such as cell type, culture condition, day captured, etc.)
#featureData, an AnnotatedDataFrame object, where rows are features (e.g. genes), and columns are gene attributes, such as biotype, gc content, etc.

#The expression value matrix must:
#have the same number of columns as the phenoData has rows.
#have the same number of rows as the featureData data frame has rows.

#Additionally:
#row names of the phenoData object should match the column names of the expression matrix.
#row names of the featureData object should match row names of the expression matrix.
#one of the columns of the featureData should be named "gene_short_name".
#this function was built for gene names, so convert protein IDs to gene Symbols from ENSEMBL

expr_matrix <- counts(spe) #exprs matrix

#convert to sparse matrix class dgMatrix
#expr_matrix <- Matrix(expr_matrix, sparse = TRUE)
#writeMM(expr_matrix, "expr_matrix.mtx")
#write.table(expr_matrix,"expr_matrix.txt",sep="\t",row.names=FALSE)

pData <- colData(spe) #cells and their metadata
write.table(pData,"pdata.txt",sep="\t",row.names=FALSE)
pdata <- read.delim("pdata.txt")
fdata <- read.delim("fdata.txt") # metadata of proteins - change column name to gene_short_name in the original panel.csv file

#pd <- new("AnnotatedDataFrame", data = pdata)
#fd <- new("AnnotatedDataFrame", data = fdata)

colnames(expr_matrix) <- row.names(pdata) #this is important
row.names(expr_matrix) <- fdata$gene_short_name
row.names(fdata) <- row.names(expr_matrix) #this is important

spe2 <- new_cell_data_set(as(expr_matrix, "dgCMatrix"),
                          cell_metadata = pdata,
                          gene_metadata = fdata)

saveRDS(spe2, "spe_CDS.rds")
spe2 <- readRDS("spe_CDS.rds")

#3.2.2- Constructing a marker file
system.file(package = "garnett")
marker_file_path <- system.file("extdata", "markers.txt", #save marker files in this directory
                                package = "garnett")


#3.2.3- Checking your markers
#Because defining the marker file is often the hardest part of the process, 
#Garnett includes functions to check whether your markers are likely to work well. 
#The two functions relevant are check_markers and plot_markers. 
#check_markers generates a table of information about your markers and plot_markers plots the most relevant information.
#In addition to the small included dataset, we have included two example marker files with the package. 

marker_check <- check_markers(cds = spe2, marker_file_path, #check markers file in the Visual Studio Code if there are errors wity typing
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")

plot_markers(marker_check)


#3.2.4- Train the classifier
#Now it's time to train the classifier. 
#The arguments should be pretty close to those for check_markers. 
#The one parameter I am changing from default below is the num_unknown argument. 
#This tells Garnett how many outgroup cells it should compare against. 
#The default is 500, but in this toy dataset with so few cells we want fewer.

set.seed(260)
classifier <- train_cell_classifier(cds = spe2,
                                    marker_file = marker_file_path,
                                    db=org.Hs.eg.db,
                                    cds_gene_id_type = "SYMBOL", #convert protein symbols to geneID symbols in the original files generating the cds object
                                    num_unknown = 50,
                                    marker_file_gene_id_type = "SYMBOL", cores = 4)

saveRDS(classifier, "Garnett_classifier.rds")

#3.2.5- Viewing the classification genes
#Garnett classification is trained using a multinomial elastic-net regression. 
#This means that certain genes are chosen as the relevant genes for distinguishing between cell types. 
#Which genes are chosen may be of interest, so Garnett includes a function to access the chosen genes. 
#Note: Garnett does not regularize the input markers, so they will be included in the classifier regardless.
#The function we use to see the relevant genes is get_feature_genes. 
#The arguments are the classifier, which node you'd like to view (if your tree is hierarchical) - use "root" for the top node and the parent cell type name for other nodes, and the db for your species.
#The function will automatically convert the gene IDs to SYMBOL if you set convert_ids = TRUE.

feature_genes <- get_feature_genes(classifier, 
                                   node = "root",
                                   db = org.Hs.eg.db, convert_ids = TRUE)
head(feature_genes)



#3.2.6- Classifying cells

spe2 <- classify_cells(spe2, classifier,
                       db = org.Hs.eg.db,
                       cluster_extend = TRUE,
                       cds_gene_id_type = "SYMBOL")

head(pData(spe2))
table(pData(spe2)$cell_type)
table(pData(spe2)$cluster_ext_type)
saveRDS(spe2, "spe_CDS.rds")

# 3.2.7-Store labels of classified cells in Spatial Experiment object

#Define color code for each cluster label
celltype <- setNames(c("#3F1B03", "#F4AD31", "#894F36", "#1C750C", "#EF8ECC", 
                       "#6471E2", "#4DB23B", "grey", "#F4800C", "#BF0A3D", "#066970", "black"),
                     c("Tumor", "Stroma", "Myeloid", "CD8 T cell", "Plasma cell", 
                       "Treg", "CD4 T cell", "Unknown", "BnT cells", "B cell", "Neutrophil", "T cell"))

metadata(spe)$color_vectors$celltype <- celltype

cell_labels <- pData(spe2)$cluster_ext_type
spe$celltype <- cell_labels
table(spe$celltype, spe$patient_id)
#check whether the cell label column was added to the cell metadata
colData(spe)


#Finally, we save the updated SpatialExperiment object.
saveRDS(spe, "spe.rds")
