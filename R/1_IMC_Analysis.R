library(imcRtools)
library(cytomapper)
library(openxlsx)
library(stringr)
library(dittoSeq)
library(RColorBrewer)

# 1- Read in single-cell information
# 1.1- steinbock generated data
spe <- read_steinbock("steinbock_outputs/")
spe

#By default, single-cell data is read in as SpatialExperiment object. 
#The summarized pixel intensities per channel and cell (here mean intensity) are stored in the counts slot. 
#Columns represent cells and rows represent channels.
counts(spe)[1:5,1:5]

#Metadata associated to individual cells are stored in the colData slot. 
#After initial image processing, these metadata include the numeric identifier (ObjectNumber), 
#the area, and morphological features of each cell. 
#In addition, sample_id stores the image name from which each cell was extracted and the width and height of the corresponding images are stored.
head(colData(spe))

#The main difference between the SpatialExperiment and the SingleCellExperiment data container in the current setting is the way spatial locations of all cells are stored. 
#For the SingleCellExperiment container, the locations are stored in the colData slot 
#while the SpatialExperiment container stores them in the spatialCoords slot:
head(spatialCoords(spe))

#The spatial object graphs generated by steinbock are read into a colPair slot of the SpatialExperiment (or SingleCellExperiment) object. 
#Cell-cell interactions (cells in close spatial proximity) are represented as “edge list” (stored as SelfHits object). 
#Here, the left side represents the column indices of the “from” cells and the right side represents the column indices of the “to” cells.
colPair(spe, "neighborhood")

#Finally, metadata regarding the channels are stored in the rowData slot. 
#This information is extracted from the panel.csv file. 
#Channels are ordered by isotope mass and therefore match the channel order of the multi-channel images
head(rowData(spe))

#1.2- IMC Segmentation Pipeline generated data
spe2 <- read_cpout("imc_segmentation_pipeline_outputs/analysis/cpout/")
rownames(spe2) <- rowData(spe2)$Clean_Target
spe2

#cell morphological features and image level metadata
head(colData(spe2))

#interaction information
colPair(spe2, type = "neighborhood")

#mean intensity per channel and cell
counts(spe2)

#2- Single-cell processing
#After reading in the single-cell data, few further processing steps need to be taken.

#2.1- Add aditional metadata
#We can set the colnames of the object to generate unique identifiers per cell:
colnames(spe) <- paste0(spe$sample_id, "_", spe$ObjectNumber)
colnames(spe)

#It is also often the case that sample-specific metadata are available externally. 
#For the current data, we need to link the cancer type (also referred to as “Indication”) to each sample. 
#This metadata is available as external excel file:
library(openxlsx)
library(stringr)
meta <- read.xlsx("steinbock_outputs/sample_metadata.xlsx")

spe$patient_id <- as.vector(str_extract_all(spe$sample_id, "Patient[1-9]", simplify = TRUE))
spe$ROI <- as.vector(str_extract_all(spe$sample_id, "00[1-3]", simplify = TRUE))
spe$indication <- meta$Indication[match(spe$patient_id, meta$Sample.ID)]

unique(spe$indication)
unique(spe$ROI)
unique(spe$patient_id)

#2.2- Transform counts - similar to quantile normalization
#The distribution of expression counts across cells is often observed to be skewed towards the right side meaning lots of cells display low counts and few cells have high counts.
#To avoid analysis biases from these high-expressing cells, the expression counts are commonly transformed or clipped.
#Here, we perform counts transformation using an inverse hyperbolic sine function. This transformation is commonly applied to flow cytometry data. 
#The cofactor here defines the expression range on which no scaling is performed. 
#While the cofactor for CyTOF data is often set to 5, IMC data usually display much lower counts. We therefore apply a cofactor of 1.
#However, other transformations such as log(counts(spe) + 0.01) should be tested when analysing IMC data.

library(dittoSeq)
dittoRidgePlot(spe, var = "CD3", group.by = "patient_id", assay = "counts") +
  ggtitle("CD3 - before transformation")

#the transformation step
assay(spe, "exprs") <- asinh(counts(spe)/1)

dittoRidgePlot(spe, var = "CD3", group.by = "patient_id", assay = "exprs") +
  ggtitle("CD3 - after transformation")

#2.3- Define interesting channels
#For downstream analysis such as visualization, dimensionality reduction and clustering, only a subset of markers should be used. 
#As convenience, we can store an additional entry in the rowData slot that specifies the markers of interest. 
#Here, we deselect the nuclear markers, which were primarily used for cell segmentation, and keep all other biological targets.

rowData(spe)$use_channel <- !grepl("DNA|Histone", rownames(spe))

#2.4- Define color schemes
#We will define color schemes for different metadata entries of the data and conveniently store them 
#in the metadata slot of the SpatialExperiment which will be helpful for downstream data visualizations. 
#We will use colors from the RColorBrewer and dittoSeq package but any other coloring package will suffice.

library(RColorBrewer)
color_vectors <- list()

ROI <- setNames(brewer.pal(length(unique(spe$ROI)), name = "BrBG"), 
                unique(spe$ROI))
patient_id <- setNames(brewer.pal(length(unique(spe$patient_id)), name = "Set1"), 
                       unique(spe$patient_id))
sample_id <- setNames(dittoColors(reps = 1)[seq_along(unique(spe$sample_id))], 
                      unique(spe$sample_id))
indication <- setNames(brewer.pal(length(unique(spe$indication)), name = "Set2"), 
                       unique(spe$indication))

color_vectors$ROI <- ROI
color_vectors$patient_id <- patient_id
color_vectors$sample_id <- sample_id
color_vectors$indication <- indication
color_vectors
metadata(spe)$color_vectors <- color_vectors

#3- Read in images
#The cytomapper package allows multi-channel image handling and visualization within the Bioconductor framework. 
#The most common data format for multi-channel images or segmentation masks is the TIFF file format, 
#which is used by steinbock and the IMC segementation pipeline to save images.
#Here, we will read in multi-channel images and segmentation masks into a CytoImageList data container. 
#It allows storing multiple multi-channel images and requires matched channels across all images within the object.
#The loadImages function is used to read in processed multi-channel images and their corresponding segmentation masks. 
#Of note, the multi-channel images generated by steinbock are saved as 32-bit images while the segmentation masks are saved as 16-bit images. 
#To correctly scale pixel values of the segmentation masks when reading them in set as.is = TRUE.

#load multi-channel TIFF images
images <- loadImages("steinbock_outputs/img/")

#load segmentation masks
masks <- loadImages("steinbock_outputs/masks_deepcell/", as.is = TRUE)

#In the case of multi-channel images, it is beneficial to set the channelNames for easy visualization. 
#Using the steinbock framework, the channel order of the single-cell data matches the channel order of the multi-channel images.
#However, it is recommended to make sure that the channel order is identical between the single-cell data and the images.
channelNames(images) <- rownames(spe)
images

all.equal(names(images), names(masks))

patient_id <- str_extract_all(names(images), "Patient[1-9]", simplify = TRUE)
indication <- meta$Indication[match(patient_id, meta$Sample.ID)] 

mcols(images) <- mcols(masks) <- DataFrame(sample_id = names(images),
                                           patient_id = patient_id,
                                           indication = indication)

#4- Generate single-cell data from images - MIGHT HAVE TO USE THIS
#An alternative way of generating a SingleCellExperiment object directly from the multi-channel images and 
#segmentation masks is supported by the measureObjects function of the cytomapper package. 
#For each cell present in the masks object, 
#the function computes the mean pixel intensity per channel as well as morphological features 
#(area, radius, major axis length, eccentricity) and the location of cells:

cytomapper_sce <- measureObjects(masks, image = images, img_id = "sample_id", return_as = "spe") #return_as	single character specifying the class of the returned object. This is either "sce" to return a SingleCellExperiment (default) or "spe" to return a SpatialExperiment object.
cytomapper_sce

#5- Save objects
#Finally, the generated data objects can be saved for further downstream processing and analysis.
saveRDS(spe, "spe.rds")
saveRDS(images, "images.rds")
saveRDS(masks, "masks.rds")


images <- readRDS("images.rds")
spe <- readRDS("steinbock_outputs/spe.rds")





















