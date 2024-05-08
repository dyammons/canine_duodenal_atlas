#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")

### Analysis note: 
# This script transfers the labels from the annotated dataset to the individual pre-integrated datasets and plots
# key DEGs identified throughout the manuscript. Generates Supplemental figure 4.

############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Run plotting for each sample ######## <<<<<<<<<<<<<<
############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
seu.obj.all <- readRDS("./output/s3/canine_duodenum_annotated.rds")
levels(seu.obj.all$majorID) <- c("T cell", "Epithelial", "Myeloid", "Plasma cell", "Cycling", "Mast", "B cell")

files <- list.files(path = "./output/s1/", pattern = "S1.rds", all.files = FALSE, full.names = TRUE)[c(2:4,6:9)]
seu.obj_list <- mapply(readRDS, files)

lapply(seu.obj_list, function(seu.obj){
    metaData <- seu.obj.all$majorID[seu.obj.all$orig.ident %in% unique(seu.obj$orig.ident)]
    names(metaData) <- substr(names(metaData), 1, 18)
    seu.obj <- AddMetaData(seu.obj, metaData, col.name = "majorID")
    seu.obj$majorID <- ifelse(is.na(seu.obj$majorID), "NA", as.character(seu.obj$majorID))
    seu.obj <- subset(seu.obj, invert = T, subset = majorID == "NA")
    seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name2")

    seu.obj$majorID <- factor(seu.obj$majorID, levels = levels(seu.obj.all$majorID))
    colz <- c("#FF89B3", "#C89504", "#00ABFD", "#B983FF", "#FA7476", "#0A9B9F", "#9B8EFF")
    names(colz) <- c("T cell", "Epithelial", "Myeloid", "Plasma cell", "Cycling", "Mast", "B cell")
    colz <- colz[match(names(colz), levels(seu.obj$majorID))]
    
    pi <- DimPlot(seu.obj, 
                  reduction = "umap", 
                  group.by = "majorID",
                  pt.size = 0.25,
                  label = F,
                  label.box = F,
                  repel = F,
                  cols = c("#FF89B3", "#C89504", "#00ABFD", "#B983FF", "#FA7476", "#0A9B9F", "#9B8EFF")
    )
    p1 <- formatUMAP(plot = pi) + theme(plot.title = element_text(size= 20)) + ggtitle(levels(seu.obj$name2))
    
#     features <- c("IL1B", "S100A12", 
#                   "SOD2", "KCNK16")
    
    features <- c("PLAUR", "PLAU", 
                  "IL7R", "ALAS1")
    p2 <- prettyFeats(seu.obj = seu.obj, nrow = 2, ncol = 2, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)

    p <- p1 + p2
    ggsave(paste0("./output/allCells/", unique(seu.obj$orig.ident), "_key_feats.png", sep = ""), width = 12, height = 6)

})
############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   End plotting for each sample ######## <<<<<<<<<<<<<<
############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

