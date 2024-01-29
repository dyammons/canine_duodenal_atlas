#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")

### Analysis note: 
# This script loads in the previously processed Seurat object (./output/s3/230816_duod_h3c4_NoIntrons_res1.3_dims40_dist0.3_neigh50_S3.rds)
# then subsets on myeloid cells and generates all figures assocaited with Figure 2

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin meyloid cell preprocessing   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#read in processed "All cells" dataset
seu.obj <- readRDS("./output/s3/230816_duod_h3c4_NoIntrons_res1.3_dims40_dist0.3_neigh50_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./colorID_cie3v4.csv", groupBy = "clusterID_2_1", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name2")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))

#subset on myeloid cells
seu.obj.sub <- subset(seu.obj,
                  subset = 
                  majorID ==  "myeloid")
table(seu.obj.sub$majorID)
table(seu.obj.sub$clusterID_2_1)
table(seu.obj.sub$orig.ident)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "./output/s2/", subName = "230829_myeloid_duod_h3c4_NoIntrons_2000", 
                     preSub = T, nfeatures = 2000, vars.to.regress = "percent.mt",  k = 25, ndims = 25
                    )

#clustree to determine clus resolution
# seu.obj <- readRDS(file = "./output/s2/230829_myeloid_duod_h3c4_NoIntrons_2000")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "230829_myeloid_duod_h3c4_NoIntrons_2000", test_dims = c("40","35", "30"), algorithm = 3, prefix = "integrated_snn_res.")

#visulize the data & evaluate
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230829_myeloid_duod_h3c4_NoIntrons_2000", final.dims = 40, final.res = 2, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = F,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#### cluster 11 looks to be T cells -- remove and re-integrate
seu.obj.sub <- subset(seu.obj, invert = T,
                      subset = 
                      clusterID_sub ==  "11")
table(seu.obj.sub$clusterID_sub)
table(seu.obj.sub$orig.ident)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "./output/s2/", subName = "230829_myeloid_duod_h3c4_NoIntrons_2500", 
                     preSub = T, nfeatures = 2500, vars.to.regress = "percent.mt",  k = 25, ndims = 25, saveRDS = F
                       )

# seu.obj <- readRDS(file = "./output/s2/230717_myeloid_duod_h3c3_NoIntrons_2500_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "230829_myeloid_duod_h3c4_NoIntrons_2500", test_dims = c("40","35", "30"), algorithm = 3, prefix = "integrated_snn_res.")

#visulize the data and proceed with analysis
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230829_myeloid_duod_h3c4_NoIntrons_2500", final.dims = 40,
                       final.res = 0.5, stashID = "clusterID_sub", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.25,
                       n.neighbors = 25, assay = "integrated", saveRDS = T,
                       features = c("PTPRC", "CD3E", "CD8A", "GZMA",
                                    "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                    "CD4", "MS4A1", "S100A12","MS4A2")
                       )

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin myeloid analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in processed data
seu.obj <- readRDS("./output/s3/230829_myeloid_duod_h3c4_NoIntrons_2500_res0.5_dims40_dist0.25_neigh25_S3.rds")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))
colz.df <- read.csv("./cellColz.csv", header = F)
colz.df <- colz.df[colz.df$V2 == "myeloid", ]
outName <- "myeloid"

#stash cell type IDs
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "Monocyte", "1" = "Neutrophil", 
                                   "2" = "Eosinophil", "3" = "Macrophage",
                                   "4" = "Neutrophil", "5" = "Eosinophil",
                                   "6" = "cDC1", "7" = "IL22RA2_DC",
                                   "8" = "Eosinophil")
                       )

seu.obj$majorID_sub <- Idents(seu.obj)
seu.obj$majorID_sub <- factor(seu.obj$majorID_sub, levels = c("Eosinophil","Neutrophil","Monocyte","Macrophage","cDC1","IL22RA2_DC"))

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


### supp data - cell type gene signatures

#generate violin plots for each cluster
# vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = "230829_myeloid_duod_h3c4_NoIntrons_2500",
#                      outDir = "./output/viln/myeloid/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
#                     )

#generate violin plots for each cell type
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_sub", numOfFeats = 24, outName = "240128_myeloid_duod_h3c4_NoIntrons_2500",
                     outDir = "./output/viln/myeloid/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )

#append with surfacome data
surface.markers <- read.csv("./surface_master.csv")[ ,c("UniProt.gene", "UniProt.description", "Surfaceome.Label", "Surfaceome.Label.Source")] %>% filter(!duplicated(UniProt.gene))
cluster.markers <- read.csv("./output/viln/myeloid/240128_myeloid_duod_h3c4_NoIntrons_2500_gene_list.csv")
write.csv(cluster.markers[ ,c(7,8,2:6)] %>% left_join(surface.markers, by = c("gene" = "UniProt.gene")),
          file = "./output/supplementalData/supplemental_data_4.csv", row.names = F)


### supp data -  export data for cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = "myeloid", outDir = "./output/cb_input/", 
               markers = "./output/viln/myeloid/231022_myeloid_duod_h3c4_NoIntrons_2500_gene_list.csv", 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", "majorID", "clusterID_sub", "name2", "majorID_sub", "cellSource", "clusterID_final"), skipEXPR = F,
               test = F,
                           feats = c("AIF1", "MS4A2", "IL18BP",
                                     "SELL", "S100A12","IL1B",
                                     "DLA-DRA", "CCL14", "C1QC",
                                     "MSR1","CSF1R","CCL3",
                                     "FLT3", "BATF3", "CADM1")
                          )


### Fig supp 3a - Create UMAP by clusterID_sub
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = TRUE
) + NoLegend()
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.8, smallAxes = T)  #, labCol = majorColors.df$labCol
ggsave(paste("./output/", outName,"/", outName, "_supp3a.png", sep = ""), width = 7, height = 7)


### Fig 2a - UMAP by clusterID_final
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_final",
              pt.size = 0.25,
              cols = colz.df$V1,
              label = T,
              label.box = T,
              shuffle = TRUE
) + NoLegend()
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.8, smallAxes = T)  #, labCol = majorColors.df$labCol
ggsave(paste0("./output/", outName, "/", outName, "_fig2a.png"), width = 7, height = 7)


### Fig 2b - create violin plots for key feats
features <- c("IL22RA2","FSCN1","LY86", "BATF3", 
              "FLT3", "CADM1",              
              "MSR1","CCL3","CSF1R",
              "DLA-DRA","CCL14", "C1QC",
              "SELL", "S100A12",
              "MS4A2", "IL18BP",
              "AIF1")

pi <- VlnPlot(object = seu.obj,
              pt.size = 0,
              same.y.lims = F,
              group.by = "majorID_sub",
              combine = T,
              cols = colz.df$V1,
              stack = T,
              fill.by = "ident",
              flip = T,
              features = features
             ) + NoLegend() + theme(axis.ticks = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.title.x = element_blank())
ggsave(paste("./output/", outName, "/", outName, "_fig2b.png", sep = ""), width = 5, height =6)


### Fig extra - Plot key feats
features <- c("AIF1",
              "MS4A2", "CD4",
              "SELL", "S100A12","IL1B",
              "DLA-DRA",
              "CCL14", "C1QC",
              "MSR1","CSF1R","CCL3",
              "FLT3", "BATF3", "CADM1")

p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol =  3, features = features, 
                 color = "black", order = F, pt.size = 0.25, title.size = 14, noLegend = T)
ggsave(paste("./output/", outName, "/", outName, "_extraFeats.png", sep = ""), width = 9, height = 15)


### Fig supp 3b - umap by sample
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
pi <- formatUMAP(pi) + labs(colour="") + theme(legend.position = "top", 
                                               legend.direction = "horizontal",
                                               legend.justification = "center",
                                               legend.title=element_text(size=12),
                                               axis.title = element_blank(),
                                               panel.border = element_blank(),
                                               plot.margin = unit(c(-7, -7, -7, -7), "pt")
                                              ) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))
ggsave(paste("./output/", outName, "/", outName, "_supp3b.png", sep = ""), width =7, height = 7)


### Fig extra: reference map using Gut Atlas data
seu.gut.myeloid <- MuDataSeurat::ReadH5AD("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/myeloid_log_counts02_v2.h5ad")
#download reference from https://www.gutcellatlas.org/#datasets

#prep the reference
reference <- seu.gut.myeloid
reference <- SCTransform(reference, verbose = FALSE)
reference <- RunPCA(reference)

#do the mapping
anchors <- FindTransferAnchors(
    reference = reference,
    query = seu.obj,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims= 1:50
)

#transfer the annotations
predictions <- TransferData(anchorset = anchors, refdata = reference$annotation,
    dims = 1:50)
seu.obj <- AddMetaData(seu.obj, metadata = predictions)
seu.obj$predicted.human <- seu.obj$predicted.id

#plot the data
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "predicted.id",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F
)
pi <- formatUMAP(plot = pi)
ggsave(paste("./output/", outName, "/",outName, "_umap_Predicted_gutAtlas.png", sep = ""), width = 10, height = 7)


### Fig extra: reference map using canine PBMC Atlas data
reference <- readRDS(file = "../../k9_PBMC_scRNA/analysis/output/s3/final_dataSet_HvO.rds")
#download reference with utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_final_dataSet_HvO.rds.gz", dest = "/pwd/to/dir/final_dataSet_HvO.rds.gz")
reference[['integrated']] <- as(object = reference[['integrated']] , Class = "SCTAssay")
DefaultAssay(reference) <- "integrated"

#do the mapping
anchors <- FindTransferAnchors(
    reference = reference,
    query = seu.obj,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims= 1:50
)

#transfer the annotations
predictions <- TransferData(anchorset = anchors, refdata = reference$celltype.l3,
    dims = 1:50)
seu.obj <- AddMetaData(seu.obj, metadata = predictions)
seu.obj$predicted.canine <- seu.obj$predicted.id

#plot the data
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "predicted.id",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F
)
pi <- formatUMAP(plot = pi)
ggsave(paste("./output/", outName, "/",outName, "_umap_Predicted_canPBMC_Atlas.png", sep = ""), width = 10, height = 7)


### Fig supp 3c - evlauate cell frequency by majorID
freqy <- freqPlots(seu.obj, method = 1, nrow= 1, groupBy = "majorID_sub", legTitle = "Cell source",refVal = "name2", showPval = T,
                   namez = unique(seu.obj$name2), 
                   colz = unique(seu.obj$colz)
                  )
ggsave(paste("./output/", outName, "/",outName, "_supp3c.png", sep = ""), width = 9, height = 3)


### Fig extra - evlauate cell frequency by clusterID_sub
freqy <- freqPlots(seu.obj, method = 1, nrow= 1, groupBy = "clusterID_sub", legTitle = "Cell source",refVal = "name2", showPval = T,
                   namez = unique(seu.obj$name2), 
                   colz = unique(seu.obj$colz)
                  )
ggsave(paste("./output/", outName, "/",outName, "_freq_by_clusID.png", sep = ""), width = 9, height = 3)


### Fig supp 3d: stacked bar graph by majorID_sub
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "cellSource", groupBy = "name2", clusters = "majorID_sub") +
scale_fill_manual(labels = unique(seu.obj$name2), 
                  values = unique(seu.obj$colz)) + theme(axis.title.y = element_blank(),
                                                         axis.title.x = element_text(size = 14),
                                                         axis.text = element_text(size = 12)) + scale_x_discrete(limits=c("Eosinophil","Neutrophil","Monocyte","Macrophage","cDC1","IL22RA2_DC"),expand = c(0, 0))
ggsave(paste("./output/", outName,"/",outName, "_supp3d.png", sep = ""), width =7, height = 5)


### Fig extra - pseudobulk DEGs btwn healthy and cie
seu.obj$allCells <- "All cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
createPB(seu.obj = seu.obj, groupBy = "allCells", comp = "cellSource", biologicalRep = "name2", lowFilter = T, dwnSam =F, min.cell = 15,
                     clusters = NULL, outDir = paste0("./output/", outName,"/pseudoBulk/") , grepTerm = "H", grepLabel = c("Healthy","CIE") #improve - fix this so it is more functional
                    )
p <- pseudoDEG(metaPWD = paste0("./output/", outName,"/pseudoBulk/allCells_deg_metaData.csv"), returnDDS = F, 
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("./output/", outName,"/pseudoBulk/"), outName = "allCells", idents.1_NAME = "CIE", idents.2_NAME = "Healthy",
          inDir = paste0("./output/", outName,"/pseudoBulk/"), title = "All cells", fromFile = T, meta = NULL, pbj = NULL, returnVolc = T, paired = F, pairBy = "", 
          minimalOuts = F, saveSigRes = T, filterTerm = "^ENSCAF", addLabs = NULL, mkDir = T
                     )


### Fig 2c - deg scatter plot
seu.obj$allCells <- "DGE analysis of myeloid cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
linDEG(seu.obj = seu.obj, groupBy = "allCells", comparision = "cellSource", outDir = paste0("./output/", outName,"/"), 
       outName = "fig2c", labCutoff = 10, contrast = c("CIE", "Healthy"),
       subtitle = T, pValCutoff = 0.01, saveGeneList = T, addLabs = ""
      )


### Fig supp 3e - gsea plot using degs
res.df <- read.csv("./output/myeloid/fig2c_DGE_analysis_of_myeloid_cells_geneList.csv")
geneList <- res.df %>% filter(p_val < 0.01) %>% filter(avg_log2FC > 0) %>% pull(X)

p <- plotGSEA(geneList = geneList, upOnly = T, category = "C2", subcategory = "CP:REACTOME", size = 4.5,termsTOplot = 10) 
pi <- p + scale_x_continuous(limits = c(-10,ceiling(max(p$data$x_axis)*1.05)), breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),name = "log10(padj)") + ggtitle("Reactome") + theme(plot.title = element_text(size = 20, hjust = 0.5),
axis.title=element_text(size = 16))
ggsave(paste("./output/", outName, "/", outName, "_supp3e.png", sep = ""), width = 10, height =7)


### Fig extra - deg between cie and healthy within each cluster
linDEG(seu.obj = seu.obj, threshold = 1, thresLine = F, groupBy = "majorID_sub", comparision = "cellSource", contrast = c("CIE", "Healthy"),
       outDir = paste0("./output/", outName,"/linDEG/"), outName = "myeloid", cluster = NULL, labCutoff = 10, noTitle = F,
       colUp = "red", colDwn = "blue", subtitle = T, returnUpList = F, returnDwnList = F, forceReturn = F, useLineThreshold = F, pValCutoff = 0.01, saveGeneList = T, addLabs = ""
                  )


### Fig 2d - deg between cie and healthy within each cluster
set.seed(12)
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(seu.obj, downsample = min(table(seu.obj$cellSource)))
features <- c("IL1A", "IL1B", "S100A12" , "SOD2", "CXCL8")
p <- FeaturePlot(seu.obj.sub,features = features, pt.size = 0.1, split.by = "cellSource", order = T, by.col = F,
                ) + labs(x = "UMAP1", y = "UMAP2") & theme(axis.text = element_blank(),
                                                           axis.title.y.right = element_text(size = 16),
                                                           axis.ticks = element_blank(),
                                                           axis.title = element_blank(),
                                                           axis.line = element_blank(),
                                                           plot.title = element_text(size=16),
                                                           title = element_blank(),
                                                           plot.margin = unit(c(0, 0, 0, 0), "cm")
                                                          ) & scale_color_gradient(breaks = pretty_breaks(n = 3), 
                                                                                   limits = c(NA, NA), low = "lightgrey", 
                                                                                   high = "darkblue")
ggsave(paste("./output/", outName, "/",outName, "_fig2d.png", sep = ""), width = 6, height = 3)


########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   emd myeloid analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<