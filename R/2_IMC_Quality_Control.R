library(imcRtools)
library(cytomapper)
library(openxlsx)
library(stringr)
library(dittoSeq)
library(RColorBrewer)

#6- Image and cell-level quality control
#The following section discusses possible quality indicators for data obtained by IMC 
#and other highly multiplexed imaging technologies. 
#Here, we will focus on describing quality metrics on the single-cell as well as image level.

#6.1- Read in the data
#We will first read in the data processed in previous sections:
images <- readRDS("images.rds")
masks <- readRDS("masks.rds")
spe <- readRDS("spe.rds")

#6.2- Segmentation quality control
library(cytomapper)
set.seed(20220118)
img_ids <- sample(seq_len(length(images)), 3)

# Normalize and clip images
cur_images <- images[img_ids]
cur_images <- normalize(cur_images, separateImages = TRUE)
cur_images <- normalize(cur_images, inputRange = c(0, 0.2))

plotPixels(cur_images,
           mask = masks[img_ids],
           img_id = "sample_id",
           missing_colour = "white",
           colour_by = c("CD163", "CD20", "CD3", "Ecad", "DNA1"),
           colour = list(CD163 = c("black", "yellow"),
                         CD20 = c("black", "red"),
                         CD3 = c("black", "green"),
                         Ecad = c("black", "cyan"),
                         DNA1 = c("black", "blue")),
           image_title = NULL,
           legend = list(colour_by.title.cex = 0.7,
                         colour_by.labels.cex = 0.7))

#An additional approach to observe cell segmentation quality and potentially also antibody specificity issues is 
#to visualize single-cell expression in form of a heatmap. 
#Here, we sub-sample the dataset to 2000 cells for visualization purposes 
#and overlay the cancer type from which the cells were extracted.
library(dittoSeq)
library(viridis)
cur_cells <- sample(seq_len(ncol(spe)), 2000)
dittoHeatmap(spe[,cur_cells], genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "exprs", cluster_cols = TRUE, scale = "none",
             heatmap.colors = viridis(100), annot.by = "indication",
             annotation_colors = list(indication = metadata(spe)$color_vectors$indication))

#6.3- Image-level quality control
#The plot below shows the
#average SNR versus the average signal intensity across all images.
library(tidyverse)
library(ggrepel)
library(EBImage)
cur_snr <- lapply(images, function(img){
  mat <- apply(img, 3, function(ch){
    # Otsu threshold
    thres <- otsu(ch, range = c(min(ch), max(ch)))
    # Signal-to-noise ratio
    snr <- mean(ch[ch > thres]) / mean(ch[ch <= thres])
    # Signal intensity
    ps <- mean(ch[ch > thres])
    
    return(c(snr = snr, ps = ps))
  })
  t(mat) %>% as.data.frame() %>% 
    mutate(marker = colnames(mat)) %>% 
    pivot_longer(cols = c(snr, ps))
})

cur_snr <- do.call(rbind, cur_snr)

cur_snr %>% 
  group_by(marker, name) %>%
  summarize(mean = mean(value),
            ci = qnorm(0.975)*sd(value)/sqrt(n())) %>%
  pivot_wider(names_from = name, values_from = c(mean, ci)) %>%
  ggplot() +
  #    geom_errorbar(aes(y = log2(mean_snr), xmin = log2(mean_ps - ci_ps), 
  #                      xmax = log2(mean_ps + ci_ps))) +
  #    geom_errorbar(aes(x = log2(mean_ps), ymin = log2(mean_snr - ci_snr), 
  #                      ymax = log2(mean_snr + ci_snr))) +
  geom_point(aes(log2(mean_ps), log2(mean_snr))) +
  geom_label_repel(aes(log2(mean_ps), log2(mean_snr), label = marker)) +
  theme_minimal(base_size = 15) + ylab("Signal-to-noise ratio [log2]") +
  xlab("Signal intensity [log2]")

#Another quality indicator is the image area covered by cells (or biological tissue). 
#This metric identifies ROIs where little cells are present, possibly hinting at incorrect selection of the ROI. 
#We can compute the percentage of covered image area using the metadata contained in the SpatialExperiment object:
colData(spe) %>%
  as.data.frame() %>%
  group_by(sample_id) %>%
  summarize(cell_area = sum(area),
            no_pixels = mean(width_px) * mean(height_px)) %>%
  mutate(covered_area = cell_area / no_pixels) %>%
  ggplot() +
  geom_point(aes(sample_id, covered_area)) + 
  theme_minimal(base_size = 15) +
  ylim(c(0, 1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ylab("% covered area") + xlab("")

#We observe that two of the 14 images show unusually low cell coverage. 
#These two images can now be visualized using cytomapper.
# Normalize and clip images
cur_images <- images[c("Patient4_005", "Patient4_007")]
cur_images <- normalize(cur_images, separateImages = TRUE)
cur_images <- normalize(cur_images, inputRange = c(0, 0.2))

plotPixels(cur_images,
           mask = masks[c("Patient4_005", "Patient4_007")],
           img_id = "sample_id",
           missing_colour = "white",
           colour_by = c("CD163", "CD20", "CD3", "Ecad", "DNA1"),
           colour = list(CD163 = c("black", "yellow"),
                         CD20 = c("black", "red"),
                         CD3 = c("black", "green"),
                         Ecad = c("black", "cyan"),
                         DNA1 = c("black", "blue")),
           legend = list(colour_by.title.cex = 0.7,
                         colour_by.labels.cex = 0.7))

#Finally, it can be beneficial to visualize the mean marker expression per image to identify images with outlying marker expression. 
#This check does not indicate image quality per se but can highlight biological differences. 
#Here, we will use the aggregateAcrossCells function of the scuttle package to compute the mean expression per image. 
#For visualization purposes, we again asinh transform the mean expression values.
library(scuttle)
image_mean <- aggregateAcrossCells(spe, 
                                   ids = spe$sample_id, 
                                   statistics="mean",
                                   use.assay.type = "counts")
assay(image_mean, "exprs") <- asinh(counts(image_mean))

dittoHeatmap(image_mean, genes = rownames(spe)[rowData(spe)$use_channel],
             assay = "exprs", cluster_cols = TRUE, scale = "none",
             heatmap.colors = viridis(100), 
             annot.by = c("indication", "patient_id", "ROI"),
             annotation_colors = list(indication = metadata(spe)$color_vectors$indication,
                                      patient_id = metadata(spe)$color_vectors$patient_id,
                                      ROI = metadata(spe)$color_vectors$ROI),
             show_colnames = TRUE)


#6.4-Cell-level quality control
#We calculate the SNR and signal intensity by fitting the mixture model across the transformed counts of all cells
#contained in the SpatialExperiment object.
library(mclust)

set.seed(220224)
mat <- apply(assay(spe, "exprs"), 1, function(x){
  cur_model <- Mclust(x, G = 2)
  mean1 <- mean(x[cur_model$classification == 1])
  mean2 <- mean(x[cur_model$classification == 2])
  
  signal <- ifelse(mean1 > mean2, mean1, mean2)
  noise <- ifelse(mean1 > mean2, mean2, mean1)
  
  return(c(snr = signal/noise, ps = signal))
})

cur_snr <- t(mat) %>% as.data.frame() %>% 
  mutate(marker = colnames(mat))

cur_snr %>% ggplot() +
  geom_point(aes(log2(ps), log2(snr))) +
  geom_label_repel(aes(log2(ps), log2(snr), label = marker)) +
  theme_minimal(base_size = 15) + ylab("Signal-to-noise ratio [log2]") +
  xlab("Signal intensity [log2]")

#Next, we observe the distributions of cell size across the individual images. 
#Differences in cell size distributions can indicate segmentation biases due to differences 
#in cell density or can indicate biological differences due to cell type compositions 
#(tumor cells tend to be larger than immune cells).
colData(spe) %>%
  as.data.frame() %>%
  group_by(sample_id) %>%
  ggplot() +
  geom_boxplot(aes(sample_id, area)) +
  theme_minimal(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ylab("Cell area") + xlab("")

summary(spe$area)

#We detect very small cells in the dataset and will remove them. 
#The chosen threshold is arbitrary and needs to be adjusted per dataset. CAREFUL HERE
sum(spe$area < 5)
spe <- spe[,spe$area >= 5]

#Another quality indicator can be an absolute measure of cell density often reported in cells per mm2
colData(spe) %>%
  as.data.frame() %>%
  group_by(sample_id) %>%
  summarize(cell_count = n(),
            no_pixels = mean(width_px) * mean(height_px)) %>%
  mutate(cells_per_mm2 = cell_count/(no_pixels/1000000)) %>%
  ggplot() +
  geom_point(aes(sample_id, cells_per_mm2)) + 
  theme_minimal(base_size = 15)  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  ylab("Cells per mm2") + xlab("")

#check differences in staining patterns
multi_dittoPlot(spe, vars = rownames(spe)[rowData(spe)$use_channel],
                group.by = "patient_id", plots = c("ridgeplot"), 
                assay = "exprs", 
                color.panel = metadata(spe)$color_vectors$patient_id)

#Finally, we will use non-linear dimensionality reduction methods to project cells from a high-dimensional (40) 
#down to a low-dimensional (2) space. 
#For this the scater package provides the runUMAP and runTSNE function. 
#To ensure reproducibility, we will need to set a seed; however different seeds and different parameter settings 
#(e.g. the perplexity) parameter in the runTSNE function need to be tested to avoid interpreting visualization artefacts. 
#For dimensionality reduction, we will use all channels that show biological variation across the dataset. 
#However, marker selection can be performed with different biological questions in mind.
library(scater)
set.seed(220225)
spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs") 
spe <- runTSNE(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs") 

#After dimensionality reduction, the low-dimensional embeddings are stored in the reducedDim slot.
reducedDims(spe)
head(reducedDim(spe, "UMAP"))

#Visualization of the low-dimensional embedding facilitates assessment of potential “batch effects”. 
#The dittoDimPlot function allows flexible visualization. 
#It returns ggplot objects which can be further modified.

library(patchwork)
# visualize patient id 
p1 <- dittoDimPlot(spe, var = "patient_id", reduction.use = "UMAP", size = 0.2) + 
  scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
  ggtitle("Patient ID on UMAP")
p2 <- dittoDimPlot(spe, var = "patient_id", reduction.use = "TSNE", size = 0.2) + 
  scale_color_manual(values = metadata(spe)$color_vectors$patient_id) +
  ggtitle("Patient ID on TSNE")

# visualize region of interest id
p3 <- dittoDimPlot(spe, var = "ROI", reduction.use = "UMAP", size = 0.2) + 
  scale_color_manual(values = metadata(spe)$color_vectors$ROI) +
  ggtitle("ROI ID on UMAP")
p4 <- dittoDimPlot(spe, var = "ROI", reduction.use = "TSNE", size = 0.2) + 
  scale_color_manual(values = metadata(spe)$color_vectors$ROI) +
  ggtitle("ROI ID on TSNE")

# visualize indication
p5 <- dittoDimPlot(spe, var = "indication", reduction.use = "UMAP", size = 0.2) + 
  scale_color_manual(values = metadata(spe)$color_vectors$indication) +
  ggtitle("Indication on UMAP")
p6 <- dittoDimPlot(spe, var = "indication", reduction.use = "TSNE", size = 0.2) + 
  scale_color_manual(values = metadata(spe)$color_vectors$indication) +
  ggtitle("Indication on TSNE")

(p1 + p2) / (p3 + p4) / (p5 + p6)

# visualize marker expression
p1 <- dittoDimPlot(spe, var = "Ecad", reduction.use = "UMAP", 
                   assay = "exprs", size = 0.2) +
  scale_color_viridis(name = "Ecad") +
  ggtitle("E-Cadherin expression on UMAP")
p2 <- dittoDimPlot(spe, var = "CD45RO", reduction.use = "UMAP", 
                   assay = "exprs", size = 0.2) +
  scale_color_viridis(name = "CD45RO") +
  ggtitle("CD45RO expression on UMAP")
p3 <- dittoDimPlot(spe, var = "Ecad", reduction.use = "TSNE", 
                   assay = "exprs", size = 0.2) +
  scale_color_viridis(name = "Ecad") +
  ggtitle("Ecad expression on TSNE")
p4 <- dittoDimPlot(spe, var = "CD45RO", reduction.use = "TSNE", 
                   assay = "exprs", size = 0.2) +
  scale_color_viridis(name = "CD45RO") +
  ggtitle("CD45RO expression on TSNE")

(p1 + p2) / (p3 + p4)

#The modified SpatialExperiment object is saved for further downstream analysis.
saveRDS(spe, "spe.rds")
