#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")

### Analysis note 1: 
# The published dataset consisted of 3 healthy and 4 CIE affected dogs. During data collection we obatined data from an
# additional 2 samples that were later determined to not have CIE. Our starting point in analysis begins with QC filtering
# using all 9 samples, but we did not end up releasing the 2 excluded samples on NCBI GEO due to being unincoperated in the
# manuscript. If you have concerns about this initial scripts and want access to the raw data, please email Dylan Ammons
# at dylan.ammons@colostate.edu and/or submit an issue on the GitHub page. The count matrices and processed data that 
# is generated from this script is provided on Zenodo. Please see the input README.md file in this repository for download instructions.

### Analysis note 2: 
# Several custom fuctions are used thoughout the analysis scripts provided. The source file for custom functions is available in customFunctions.R

### Analysis note 3: 
# This script loads in 10x count matrices, filters them linently then increases stringency on non-epithieal cells.
# The output is a Seurat object saved as ./output/s3/230706_duod_h3c6_NoIntrons_res0.4_dims50_dist0.5_neigh40_S3.rds

# Thank you and please cite our work if you use our data in your project. And feel free to reach out with questions/concerns.

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin all cell preprocessing   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in 10x data and qc filter eeach sample
load10x(din = "./inputNoIntrons/", dout = "./output/s1/", outName = "230403_NoIntrons", testQC = F, nFeature_RNA_high = 3000, nFeature_RNA_low = 200, percent.mt_high = 80, nCount_RNA_high = 25000, nCount_RNA_low = 100)

#integrate the data
seu.obj <- sctIntegrate(din = "./output/s1/", dout = "./output/s2/", outName = "230706_duod_h3c6_NoIntrons", vars.to.regress = "percent.mt", nfeatures = 2500)

#plot clustree to idetify cluster res to use
seu.obj <- readRDS("./output/s2/230706_duod_h3c6_NoIntrons_seu.integrated.obj_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "230706_duod_h3c6_NoIntrons", test_dims = c(50), algorithm = 3, prefix = "integrated_snn_res.")

#plot the data
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230706_duod_h3c6_NoIntrons", final.dims = 50, final.res = 0.4, stashID = "clusterID_pre", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.4, n.neighbors = 30, assay = "integrated", saveRDS = F, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#stash addtional metadata regarding diease state
Idents(seu.obj) <- "orig.ident"
seu.obj$cellSource <- ifelse(grepl("norm", seu.obj@meta.data$orig.ident), "Healthy", "CIE")

#check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", outName, "_QC_feats.png", sep = ""), width = 9, height = 3)

#ID non-duodenal cells to complete secondary pct.mt filter
duod <- c(2,3,5,14,17,10,21,18)

#remove non-duod cells with pct.mt > 12.5
seu.obj$quality <- ifelse(!seu.obj$clusterID_pre %in% duod & seu.obj$percent.mt > 12.5, "lowQuality", "fairQuality")
seu.obj <- subset(seu.obj,
                  subset = quality != "lowQuality"
                 )


#plot clustree to idetify cluster res to use
DefaultAssay(seu.obj) <- "integrated"
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "230706_duod_h3c6_NoIntrons", test_dims = c(50), algorithm = 3, prefix = "integrated_snn_res.")

#repeat clustering after additonal filtering
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230706_duod_h3c6_NoIntrons", final.dims = 50, final.res = 0.4, stashID = "clusterID", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.5, n.neighbors = 40, assay = "integrated", saveRDS = F, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#save the RDS for downstream analysis
saveRDS(seu.obj, file = "./output/s3/230706_duod_h3c6_NoIntrons_res0.4_dims50_dist0.5_neigh40_S3.rds")

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end all cell preprocessing   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<