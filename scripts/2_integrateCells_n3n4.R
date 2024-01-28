#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")

### Analysis note: 
# This script loads in the previously processed Seurat object (./output/s3/230706_duod_h3c6_NoIntrons_res0.4_dims50_dist0.5_neigh40_S3.rds)
# then removes the two dogs who ended up not meeting inclusion criteria. The data are re-integrated then 1 low quality cluster is removed.
# Output is a processed Seurat object called ./output/s3/230816_duod_h3c4_NoIntrons_res1.3_dims40_dist0.3_neigh50_S3.rds

################################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin all cell preprocessing   ######## <<<<<<<<<<<<<<
################################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in the object that contains non-CIE diseased dogs
seu.obj <- readRDS("./output/s3/230706_duod_h3c6_NoIntrons_res0.4_dims50_dist0.5_neigh40_S3.rds")

#load in meta data
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "exclude")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy", "CIE"))
outName <- "allCells"
colArray <- read.csv("./refColz.csv", header = T)

#subset on samples to include
table(seu.obj$exclude)
seu.obj <- subset(seu.obj,
                  subset = 
                  exclude ==  "keep")

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = "230706_duod_h3c4_NoIntrons", preSub = T, nfeatures = 2500,
                      vars.to.regress = "percent.mt"
                       )

#clustree to determine clus resolution
# seu.obj <- readRDS(file = "./output/s2/230706_duod_h3c4_NoIntrons_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "230706_duod_h3c4_NoIntrons", test_dims = 40, algorithm = 3, prefix = "integrated_snn_res.")

#visulize the data & evaluate
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230706_duod_h3c4_NoIntrons", final.dims = 40, final.res = 0.3, stashID = "clusterID2", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 50, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


### Note: cluster 10 looks susspect -- filtering out (low nFeature and no clear defining feats)
seu.obj <- subset(seu.obj, invert = T,
                  subset = 
                  clusterID2 ==  "10")
table(seu.obj$clusterID2)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = "230816_duod_h3c4_NoIntrons", preSub = T, nfeatures = 2500,
                      vars.to.regress = "percent.mt"
                       )

#clustree to determine clus resolution
# seu.obj <- readRDS(file = "./output/s2/230713_duod_h3c4_NoIntrons_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "230606_duod_h3c3_NoIntrons", test_dims = 40, algorithm = 3, prefix = "integrated_snn_res.")

#visulize the data & evaluate
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230816_duod_h3c4_NoIntrons", final.dims = 40, final.res = 1.3, stashID = "clusterID_2_1", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 50, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end all cell preprocessing   ######## <<<<<<<<<<<<<<
############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<