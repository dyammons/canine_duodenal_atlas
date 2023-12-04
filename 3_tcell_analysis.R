#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin T cell preprocessing   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#read in processed "All cells" dataset
seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/s3/230816_duod_h3c4_NoIntrons_res1.3_dims40_dist0.3_neigh50_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./colorID_cie3v4.csv", groupBy = "clusterID_2_1", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name2")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))
outName <- "tcell"

#subset on tcell cells
seu.obj.sub <- subset(seu.obj,
                  subset = 
                  majorID ==  "tcell")

table(seu.obj.sub$majorID)
table(seu.obj.sub$clusterID_2_1)
table(seu.obj.sub$orig.ident)


#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "./output/s2/", subName = "230829_tcell_duod_h3c4_NoIntrons_2500", 
                     preSub = T, nfeatures = 2500, vars.to.regress = "percent.mt", saveRDS = F
                    )

#clustree to determine clus resolution
# seu.obj <- readRDS(file = "./output/s2/230829_myeloid_duod_h3c4_NoIntrons_2000")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "230829_tcell_duod_h3c4_NoIntrons_2500", test_dims = c("40","35", "30"), algorithm = 3, prefix = "integrated_snn_res.")

#visulize the data & evaluate
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230829_tcell_duod_h3c4_NoIntrons_2500", final.dims = 40, final.res = 0.2, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


#generate violin plots for each cluster
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = "230829_tcell_duod_h3c4_NoIntrons_2500",
                     outDir = "./output/viln/tcell/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )


#remove susspected low quality cluster
seu.obj.sub <- subset(seu.obj, invert = T,
                  subset = 
                  clusterID_sub ==  "6")

table(seu.obj.sub$clusterID_sub)
table(seu.obj.sub$clusterID_2_1)
table(seu.obj.sub$orig.ident)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "./output/s2/", subName = "230913_tcell_duod_h3c4_NoIntrons_2500", 
                     preSub = T, nfeatures = 2500, vars.to.regress = "percent.mt", saveRDS = F
                    )

#clustree to determine clus resolution
# seu.obj <- readRDS(file = "./output/s2/230913_tcell_duod_h3c4_NoIntrons_2500_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "230913_tcell_duod_h3c4_NoIntrons_2500", test_dims = c("40","35", "30"), algorithm = 3, prefix = "integrated_snn_res.")

#visulize the data & evaluate
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230913_tcell_duod_h3c4_NoIntrons_2500", final.dims = 35, final.res = 0.6, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin T cell analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/s3/230913_tcell_duod_h3c4_NoIntrons_2500_res0.6_dims35_dist0.3_neigh30_S3.rds")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))
colz.df <- read.csv("./cellColz.csv", header = F)
outName <- "tcell"

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "Trm (c0)", "1" = "Trm (c0)", 
                                   "2" = "Tinf (c1)", "3" = "Trm (c0)",
                                   "4" = "gdT_1 (c2)", "5" = "Tinf (c1)",
                                   "6" = "NK_T (c3)", "7" = "Tinf (c1)",
                                   "8" = "CD8mem (c4)", "9" = "gdT_2 (c5)",
                                   "10" = "Trm (c0)", "11" = "Treg (c6)",
                                   "12" = "Trm (c0)", "13" = "ILC2 (c7)",
                                   "14" = "NK (c8)","15" = "T_IFN (c9)")
                       )


seu.obj$majorID_sub <- Idents(seu.obj)
seu.obj$majorID_sub <- factor(seu.obj$majorID_sub, levels = levels(seu.obj$majorID_sub)[c(1,3,6,5,4,2,8,7,9,10)])

#stash the numerical ID
clusterID_final <- table(seu.obj$majorID_sub) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
mutate(clusterID_final=row_number()-1) %>% arrange(clusterID_final) 

newID <- clusterID_final$clusterID_final
names(newID) <- clusterID_final$Var1
Idents(seu.obj) <- "majorID_sub"
seu.obj <- RenameIdents(seu.obj, newID)
table(Idents(seu.obj))
seu.obj$clusterID_final <- Idents(seu.obj)


### Fig extra - check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/",outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


### Fig extra - run singleR
singleR(seu.obj = seu.obj, outName = "230913_tcell_duod_h3c4_NoIntrons_2500", clusters = "clusterID_sub", outDir = "./output/singleR/")

#generate violin plots for each cluster
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = "230913_tcell_duod_h3c4_NoIntrons_2500",
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )

saveName <- "231022_tcell_duod_h3c4_NoIntrons_2500"
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_sub", numOfFeats = 24, outName = saveName,
                     outDir = paste0("/output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )

#export data for cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "./output/cb_input/", 
               markers = paste0("./output/viln/",outName,"/",saveName,"_gene_list.csv"), 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", "majorID",
                                                   "clusterID_sub", "name2", "majorID_sub", "cellSource", "clusterID_final"), 
               skipEXPR = F,test = F,
               feats = c("CD4", "CD8A", "GZMA", "GZMB", 
                         "CD40LG", "TCF7", "FASLG",
                         "THY1", "BIN1", "TNFRSF6B", "CTLA4",
                         "IL13", "IL17RB", "NCR3","F2RL3")
                          
                          )


### Fig Supp 2a - tcell unsupervised clustering
pi <- DimPlot(seu.obj, 
        reduction = "umap", 
        group.by = "clusterID_sub",
        pt.size = 0.25,
              cols = colz.df$V1,

        label = TRUE,
        label.box = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8) + NoLegend() + theme(axis.title = element_blank(),
                                                                                  panel.border = element_blank())

axes <- ggplot() + labs(x = "UMAP1", y = "UMAP2") + 
theme(axis.line = element_line(colour = "black", 
                               arrow = arrow(angle = 30, length = unit(0.1, "inches"),
                                             ends = "last", type = "closed"),
                              ),
      axis.title.y = element_text(colour = "black", size = 20),
      axis.title.x = element_text(colour = "black", size = 20),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
     )

p <- p + inset_element(axes,left= 0,
                       bottom = 0,
                       right = 0.25,
                       top = 0.25,
                       align_to = "full")
ggsave(paste0("./output/", outName, "/", outName, "_supp2a.png"), width = 7, height = 7)


### Fig 2a - tcell unsupervised clustering
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_final",
              pt.size = 0.25,
              cols = colz.df$V1,
              label = T,
              label.box = T
             )
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.8, labCol = c(rep("black",5),"white","black","white","white","white")) + NoLegend() + theme(axis.title = element_blank(),
      panel.border = element_blank())


axes <- ggplot() + labs(x = "UMAP1", y = "UMAP2") + 
theme(axis.line = element_line(colour = "black", 
                               arrow = arrow(angle = 30, length = unit(0.1, "inches"),
                                             ends = "last", type = "closed"),
                              ),
      axis.title.y = element_text(colour = "black", size = 12),
      axis.title.x = element_text(colour = "black", size = 12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
     )

p <- p + inset_element(axes,left= 0,
                       bottom = 0,
                       right = 0.25,
                       top = 0.25,
                       align_to = "full")
ggsave(paste0("./output/", outName, "/", outName, "_fig2a.png"), width = 7, height = 7)


### Fig extra - plot original cluster umap
pi <- DimPlot(seu.obj, 
        reduction = "umap", 
        group.by = "clusterID_2_1",
        pt.size = 0.25,
        label = TRUE,
        label.box = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8) + NoLegend()
ggsave(paste("./output/", outName,"/",outName, "_raw_orig_UMAP.png", sep = ""), width = 7, height = 7)


### Fig 2b - dot plot of key t cell features
features <- c("TRAT1","GZMA","GZMB","SYNE1","AKAP5","SEMA7A",
              "ENSCAFG00000041010","ENSCAFG00000011098","ENSCAFG00000030206","GZMK",
              "TCF7","CRTAM",
              "IL7R","VIM","CXCR4","PLAC8",
              "MMP9","IL17RB", "IL13", 
              "CTLA4","FOXP3","TNFRSF4","TNFRSF18",
              "NCR3","F2RL3","KLRB1","IL12RB2",
              "IFI44","OAS1","ISG15")


pi <- majorDot(seu.obj = seu.obj, groupBy = "majorID_sub",
               features = features
              ) + theme(axis.title = element_blank(),
                        axis.text = element_text(size = 12))
ggsave(paste("./output/", outName, "/", outName, "_fig2b.png", sep = ""), width =9, height = 5)


### Fig extra - plot on featPlots
p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 6, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 8, noLegend = T)
ggsave(paste("./output/", outName, "/",outName, "_tcell_featPlots.png", sep = ""), width = 12, height = 10)


### Fig supp 2c - transfer TRDC data from ROS data
seu.obj.ros <- readRDS(file = "./output/s3/230828_duod_h3_ros_res0.4_dims50_dist0.4_neigh30_S3.rds")

seu.obj.Hsub <- subset(seu.obj,
                  subset = 
                  cellSource ==  "Healthy"
                 ) 

colnames(seu.obj.ros@assays$RNA@data) <- ifelse(grepl("_1",rownames(seu.obj.ros@meta.data)),paste0(substr(rownames(seu.obj.ros@meta.data), 1, nchar(rownames(seu.obj.ros@meta.data))-2),"_7") ,rownames(seu.obj.ros@meta.data))
colnames(seu.obj.ros@assays$RNA@data) <- ifelse(grepl("_2",rownames(seu.obj.ros@meta.data)),paste0(substr(rownames(seu.obj.ros@meta.data), 1, nchar(rownames(seu.obj.ros@meta.data))-2),"_8") ,rownames(seu.obj.ros@meta.data))
colnames(seu.obj.ros@assays$RNA@data) <- ifelse(grepl("_3",rownames(seu.obj.ros@meta.data)),paste0(substr(rownames(seu.obj.ros@meta.data), 1, nchar(rownames(seu.obj.ros@meta.data))-2),"_9") ,rownames(seu.obj.ros@meta.data))

seu.obj.Hsub <- AddMetaData(seu.obj.Hsub, metadata = seu.obj.ros@assays$RNA@data[rownames(seu.obj.ros@assays$RNA@data) == "TRDC",], col.name = "TRDC")

features = c("TRDC", "ENSCAFG00000030206")

p <- prettyFeats(seu.obj = seu.obj.Hsub, nrow = 1, ncol = 2, features = features, color = "black", order = F) 
ggsave(paste("./output/", outName, "/", outName, "_supp_2b.png", sep = ""), width = 8, height = 4)


### Fig - reference map using PBMC data
reference <- readRDS(file = "../../k9_PBMC_scRNA/analysis/output/s3/final_dataSet_HvO.rds")
reference[['integrated']] <- as(object = reference[['integrated']] , Class = "SCTAssay")

DefaultAssay(reference) <- "integrated"

anchors <- FindTransferAnchors(
    reference = reference,
    query = seu.obj,
    normalization.method = "SCT",
    reference.reduction = "pca", #reference.reduction = "umap",
    dims= 1:50 #dims= 1:2
)

predictions <- TransferData(anchorset = anchors, refdata = reference$celltype.l3,
    dims = 1:50)
seu.obj <- AddMetaData(seu.obj, metadata = predictions)
seu.obj$predicted.canine <- seu.obj$predicted.id

pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "predicted.id",
              #cols = levels(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F
)
pi <- formatUMAP(plot = pi)
ggsave(paste("./output/", outName,"/",outName, "_umap_Predicted.png", sep = ""), width = 10, height = 7)


### Fig - reference map using gut atlas data
reference <- MuDataSeurat::ReadH5AD("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/Tcell_log_counts02_v2.h5ad")

reference <- SCTransform(reference, verbose = FALSE)
reference <- RunPCA(reference)

anchors <- FindTransferAnchors(
    reference = reference,
    query = seu.obj,
    normalization.method = "SCT",
    reference.reduction = "pca", #reference.reduction = "umap",
    dims= 1:50 #dims= 1:2
)

predictions <- TransferData(anchorset = anchors, refdata = reference$annotation,
    dims = 1:50)
seu.obj <- AddMetaData(seu.obj, metadata = predictions)
seu.obj$predicted.human <- seu.obj$predicted.id

pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "predicted.id",
              #cols = levels(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F
)
pi <- formatUMAP(plot = pi)
ggsave(paste("./output/", outName,"/",outName, "_umap_Predicted_gutAtlas.png", sep = ""), width = 10, height = 7)


### Fig supp 2c - umap by sample
Idents(seu.obj) <- "cellSource"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data$cellSource)))
pi <- DimPlot(seu.obj.ds, 
              reduction = "umap", 
              group.by = "name2",
              cols = unique(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.5,
              label = FALSE,
              shuffle = TRUE
)
pi <- formatUMAP(pi) + labs(colour="") + theme(legend.position = "top", legend.direction = "horizontal",legend.title=element_text(size=12)) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))
ggsave(paste("./output/", outName, "/", outName, "_supp2c.png", sep = ""), width =7, height = 7)


### Fig 2c - Evlauate cell frequency by majorID_sub
freqy <- freqPlots(seu.obj, method = 1, nrow= 1, groupBy = "majorID_sub", legTitle = "Cell source",refVal = "name2", showPval = T,
                   namez = unique(seu.obj$name2), 
                   colz = unique(seu.obj$colz)
                  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + NoLegend()

ggsave(paste("./output/", outName, "/",outName, "_fig2c.png", sep = ""), width = 13.5, height = 3)

#confirm variability is consistent between CIE and H
res.ftest <- lapply(levels(freqy$data$majorID_sub), function (x){
    data.df <- freqy$data[freqy$data$majorID_sub ==  x, ]
    car::leveneTest(freq ~ cellSource, data = data.df)
})


### Fig 2d - DEG between Trm and Tinf
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_final", idents.1 = "0", idents.2 = "1", bioRep = "name2",padj_cutoff = 0.05, lfcCut = 0.58, 
                      minCells = 25, outDir = paste0("./output/", outName, "/"), 
                      title = "Trm_vs_Tinf", idents.1_NAME = "Trm", idents.2_NAME = "Tinf", 
                      returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 24
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in Trm", leftLab = "Up in Tinf", arrowz = T) + labs(x = "log2(FC) Trm vs Tinf") + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_fig2d.png", sep = ""), width = 7, height = 7)


### Fig 2e: gsea of the DGE results
p <- plotGSEA(pwdTOgeneList = "./output/tcell/Trm_vs_Tinf_all_genes.csv", category = "C5", subcategory = NULL, 
              upCol = "blue", dwnCol = "red", size = 4.5) 

minVal <- -12
maxVal <- 8

pi <- p + scale_x_continuous(limits = c(minVal, maxVal), name = "Signed log10(padj)") + 
    theme(axis.title=element_text(size = 16)) + 

    annotate("segment", x = -0.1, 
             y = 17, 
             xend = minVal, 
             yend = 17, 
             lineend = "round", linejoin = "bevel", linetype ="solid", colour = "blue",
             size = 1, arrow = arrow(length = unit(0.1, "inches"))
            ) + 

    annotate(geom = "text", x = (minVal-0.1*1.5)/2-0.1*1.5, 
             y = 18,
             label = "Up in Tinf",
             hjust = 0.5,
             vjust = 1.5,
             size = 5) +

    annotate("segment", x = 0.1, 
             y = 17, 
             xend = maxVal,
             yend = 17,
             lineend = "round", linejoin = "bevel", linetype ="solid", colour = "red",
             size = 1, arrow = arrow(length = unit(0.1, "inches"))
            ) + 

    annotate(geom = "text", x = (maxVal-0.1*1.5)/2+0.1*1.5, 
             y = 18,
             label = "Up in Trm",
             hjust = 0.5,
             vjust = 1.5,
             size = 5)


ggsave(paste("./output/", outName, "/", outName, "_fig2e.png", sep = ""), width = 10, height = 7)

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end T cell analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<