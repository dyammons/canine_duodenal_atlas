#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin epithelial preprocessing   ######## <<<<<<<<<<<<<<
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

#subset on epithelial cells
seu.obj.sub <- subset(seu.obj,
                  subset = 
                  majorID ==  "epithelial" | majorID ==  "fibroblast")

table(seu.obj.sub$majorID)
table(seu.obj.sub$clusterID_2_1)
table(seu.obj.sub$orig.ident)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "./output/s2/", subName = "230916_duod_duod_h3c4_NoIntrons_2500", 
                     preSub = T, nfeatures = 2500, vars.to.regress = "percent.mt",  saveRDS = F
                    )

#clustree to determine clus resolution
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "230916_duod_duod_h3c4_NoIntrons_2500", test_dims = c("40","35", "30"), algorithm = 3, prefix = "integrated_snn_res.")

#visulize the data & evaluate
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230916_duod_duod_h3c4_NoIntrons_2500", final.dims = 40, final.res = 0.4, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = F,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#### cluster 10 looks to be T cells -- remove and re-integrate
seu.obj.sub <- subset(seu.obj, invert = T,
                      subset = 
                      clusterID_sub ==  "10")
table(seu.obj.sub$clusterID_sub)
table(seu.obj.sub$orig.ident)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "./output/s2/", subName = "230916_duod_duod_h3c4_NoIntrons_2500", 
                     preSub = T, nfeatures = 2500, vars.to.regress = "percent.mt",  saveRDS = F
                    )

clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "230916_duod_duod_h3c4_NoIntrons_2500", test_dims = c("40","35", "30"), algorithm = 3, prefix = "integrated_snn_res.")

#visulize the data and proceed with analysis
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230916_duod_duod_h3c4_NoIntrons_2500", final.dims = 40,
                       final.res = 1.7, stashID = "clusterID_sub", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.25,
                       n.neighbors = 25, assay = "integrated", saveRDS = T,
                       features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin epithelial analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in preprocessed data
seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/s3/230916_duod_duod_h3c4_NoIntrons_2500_res1.7_dims40_dist0.25_neigh25_S3.rds")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))
colz.df <- read.csv("./cellColz.csv", header = F)
colz.df <- colz.df[colz.df$V2 == "duod", ]
outName <- "duod"

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "Enterocyte_1", "1" = "Enterocyte_1", 
                                   "2" = "Enterocyte_2", "3" = "Enterocyte_1",
                                   "4" = "Enterocyte_1", "5" = "Enterocyte_2",
                                   "6" = "Enterocyte_1", "7" = "Enterocyte_1",
                                   "8" = "Enterocyte_2","9" = "Enterocyte_2",
                                  "10" = "Enterocyte_2", "11" = "Enterocyte_3", 
                                   "12" = "Enterocyte_3", "13" = "Enterocyte_1",
                                   "14" = "Enterocyte_3", "15" = "Goblet",
                                   "16" = "Tuft", "17" = "BEST4_epithelial",
                                   "18" = "IFN_enterocyte","19" = "Enterocyte_1",
                                  "20" = "Stromal", "21" = "BEST4_epithelial", 
                                   "22" = "Enterocyte_2", "23" = "Goblet",
                                   "24" = "Enteroendocrine")
                       )


seu.obj$majorID_sub <- Idents(seu.obj)
seu.obj$majorID_sub <- factor(seu.obj$majorID_sub, levels = levels(seu.obj$majorID_sub)[c(1,2,3,7,6,8,4,5,9)])

#stash the numerical ID
clusterID_final <- table(seu.obj$majorID_sub) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
mutate(clusterID_final=row_number()-1) %>% arrange(clusterID_final) 

newID <- clusterID_final$clusterID_final
names(newID) <- clusterID_final$Var1
Idents(seu.obj) <- "majorID_sub"
seu.obj <- RenameIdents(seu.obj, newID)
table(Idents(seu.obj))
seu.obj$clusterID_final <- Idents(seu.obj)

#generate violin plots for each cluster
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = "230916_duod_duod_h3c4_NoIntrons_2500",
                     outDir = "./output/viln/duod/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )

saveName <- "231022_duod_duod_h3c4_NoIntrons_2500"
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_sub", numOfFeats = 24, outName = saveName,
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )

#export data fpr cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "./output/cb_input/", 
               markers = paste0("./output/viln/",outName,"/",saveName,"_gene_list.csv"), 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", "majorID",
                                                   "clusterID_sub", "name2", "majorID_sub", "cellSource", "clusterID_final"), 
               skipEXPR = T,test = F,
               feats = c("SLC15A1", "ACE2", "IDO1", "CFTR", 
                         "BEST4", "NOTCH2", "AREG", "FABP3", 
                         "IRAG2", "SYNE2", "ONECUT2","SI")
                          
                          )


### Fig supp 4a - UMAP by clusterID_sub
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = TRUE
)
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8) + NoLegend() #, labCol = majorColors.df$labCol
ggsave(paste("./output/", outName,"/", outName, "_supp4a.png", sep = ""), width = 7, height = 7)


### Fig 4a - UMAP by majorID_sub
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_final",
              pt.size = 0.25,
              cols = colz.df$V1,
              label = T,
              label.box = T,
              shuffle = TRUE
)
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.8) + NoLegend() + theme(axis.title = element_blank(),
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
ggsave(paste0("./output/", outName, "/", outName, "_fig4a.png"), width = 7, height = 7)


### Fig 4b - dotplot of major featrures
features <- c("SLC15A1", "SI", "ACE2", "TMEM37",
              #"MGAM","SLC7A8","PTPRH","CNDP1","APOC3","MICAL3","MGST3","TPT1","COX4I1",
              "ISG15","HERC6",
              "IDO1", "CFTR", "PLCB4","BEST4","NOTCH2","ADGRF5","CEACAM20",
              "DCN","COL1A1","CALD1",
              "AREG", "SYTL2", "AGR2", "FABP3", "FOXA3", "NRTN",
              "PTGS1", "IRAG2", "FYB1", "ALOX5", "TSNARE1", "SYNE2",
              "SCG2", "SCG3", "SCG5", "RFX6", "PCSK1", "RIMBP2")

p <- majorDot(seu.obj = seu.obj, groupBy = "majorID_sub",
              yAxis = levels(seu.obj$majorID_sub),
                  features = features
                 ) + theme(axis.title = element_blank(),
                           axis.text = element_text(size = 12))
ggsave(paste("./output/", outName, "/", outName, "_fig4b.png", sep = ""), width = 9, height = 4)


### Fig supp 4b - reference map using the human gut atlas -- epithieal reference
seu.gut.duod <- MuDataSeurat::ReadH5AD("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/epi_log_counts02_v2.h5ad")
reference <- seu.gut.duod

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
pi <- formatUMAP(plot = pi) + NoLegend() + theme(plot.title = element_text(size = 18, vjust = 1),
                                                 axis.title = element_blank(),
                                                 panel.border = element_blank()) + ggtitle("Human epithelial reference mapping")
ggsave(paste("./output/", outName,"/",outName, "_umap_Predicted_gutAtlas_epi.png", sep = ""), width = 7, height = 7)


### Fig supp 4b - reference map using the human gut atlas -- mesnechymal reference
reference <- MuDataSeurat::ReadH5AD("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/Mesenchyme_log_counts02_v2.h5ad")

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
seu.obj$predicted.human.Mesenchyme <- seu.obj$predicted.id

pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "predicted.id",
              #cols = levels(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F
)
pi <- formatUMAP(plot = pi) + NoLegend() + theme(plot.title = element_text(size = 18, vjust = 1),
                                                 axis.title = element_blank(),
                                                 panel.border = element_blank()) + ggtitle("Human mesenchyme reference mapping")
ggsave(paste("./output/", outName,"/",outName, "_umap_Predicted_gutAtlas_mes.png", sep = ""), width = 7, height = 7)


### Fig supp: plot enrichment scores
ecLists <- read.csv("gut_ecTerms.csv", header = T)
ecLists <- ecLists[ecLists$lineage == "Epithelial" | ecLists$lineage == "Mesenchymal", ]

modulez <- split(ecLists$genes, ecLists$cluster)

modulez <- modulez[unname(unlist(lapply(unlist(lapply(modulez, length)), function(x){ifelse(x >= 10, TRUE, FALSE)})))]

names(modulez) <- paste0(names(modulez),"_SIG")

seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")

names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

features <- names(modulez)
features <- features[c(9,2,39,28:32,1,11,34,7,8)]

p <- majorDot(seu.obj = seu.obj, groupBy = "majorID_sub",
                     features = features
                    ) + theme(legend.position = "bottom",
                              axis.title.y = element_blank(),
                              plot.margin = margin(7, 7, 0, 100, "pt")) + scale_y_discrete(position = "right") + guides(size = guide_legend(nrow = 2, byrow = F, title = 'Percent\nenriched')) + guides(color = guide_colorbar(title = 'Scaled\nenrichment score')) 
ggsave(paste("./output/", outName, "/", outName, "_gut_modScores.png", sep = ""), width = 8, height = 6)


### Fig 4c - evlauate cell frequency by cluster
freqy <- freqPlots(seu.obj, method = 1, nrow= 2, groupBy = "majorID_sub", legTitle = "Cell source",refVal = "name2", showPval = T,
                   namez = unique(seu.obj$name2), 
                   colz = unique(seu.obj$colz)
                  ) + theme(strip.text = element_text(size=8))
ggsave(paste("./output/", outName, "/",outName, "_fig4c.png", sep = ""), width = 8.5, height = 5)

#ensure appropriate stats were run - c3 (enterocyte 3 should use a different statisical appraoch, still n.s. either way though)
res.ftest <- lapply(levels(freqy$data$majorID_sub), function (x){
    data.df <- freqy$data[freqy$data$majorID_sub ==  x, ]
    car::leveneTest(freq ~ cellSource, data = data.df)
})
res.ftest


### Fig supp 4b - pseudobulk dge analysis
seu.obj$allCells <- "Epithelial cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
createPB(seu.obj = seu.obj, groupBy = "allCells", comp = "cellSource", biologicalRep = "name2", lowFilter = T, dwnSam =F, min.cell = 15,
                     clusters = NULL, outDir = paste0("./output/", outName,"/pseudoBulk/") , grepTerm = "H", grepLabel = c("Healthy","CIE") #improve - fix this so it is more functional
                    )

p <- pseudoDEG(metaPWD = paste0("./output/", outName,"/pseudoBulk/allCells_deg_metaData.csv"), returnDDS = F, 
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("./output/", outName,"/pseudoBulk/"), outName = "allCells", idents.1_NAME = "CIE", idents.2_NAME = "Healthy",
          inDir = paste0("./output/", outName,"/pseudoBulk/"), title = "All cells", fromFile = T, meta = NULL, pbj = NULL, returnVolc = T, paired = F, pairBy = "", 
          minimalOuts = F, saveSigRes = T, filterTerm = "^ENSCAF", addLabs = NULL, mkDir = T
                     )

pi  <- prettyVolc(plot = p[[1]], rightLab = NULL, leftLab = NULL, arrowz = F) + labs(x = "log2(FC) CIE vs Normal") + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_supp4c.png", sep = ""), width = 7, height = 7)


### Fig 4d - dge scatter plot
seu.obj$allCells <- "DGE analysis of epithelial cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
linDEG(seu.obj = seu.obj, groupBy = "allCells", comparision = "cellSource", outDir = paste0("./output/", outName,"/fig4d_"), 
       outName = "epithelial", labCutoff = 10, contrast = c("CIE", "Healthy"),
       subtitle = T, pValCutoff = 0.01, saveGeneList = T, addLabs = "KCNK16"
      )

#get number of degs for reporting
df.res <- read.csv("./output/duod/epithelial_DGE analysis of epithelial cellsgeneList.csv")
df.res %>% filter(avg_log2FC > 0, p_val < 0.01) %>% nrow()
df.res %>% filter(avg_log2FC < 0, p_val < 0.01) %>% nrow()


### Fig 4e - plot KCNK16 split by CIE vs H
set.seed(12)
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(seu.obj, downsample = min(table(seu.obj$cellSource)))

features <- c("KCNK16")

p <- FeaturePlot(seu.obj.sub,features = features, pt.size = 0.01, split.by = "cellSource", order = T, by.col = F) + labs(x = "UMAP1", y = "UMAP2") & theme(axis.text = element_blank(),
                                                                                                                                                axis.title.y.right = element_text(size = 16),
                                                                                                                                            
                                                                                                                           axis.ticks = element_blank(),
                                                                                                                           axis.title = element_blank(),
                                                                                                                           axis.line = element_blank(),
                                                                                                                           plot.title = element_text(size=16),
                                                                                                                           title = element_blank(),
                                                                                                                           plot.margin = unit(c(1, 0, 0, 0), "pt")
                                                                                                                          ) & scale_color_gradient(breaks = pretty_breaks(n = 3), limits = c(NA, NA), low = "lightgrey", high = "darkblue")
ggsave(paste("./output/", outName, "/",outName, "_fig4e.png", sep = ""), width = 2, height = 4)


### Fig extra - dge scatter plot within each epithelial cell type
linDEG(seu.obj = seu.obj, groupBy = "majorID_sub", comparision = "cellSource", outDir = paste0("./output/", outName,"/"), 
       outName = "epithelial", labCutoff = 10, contrast = c("CIE", "Healthy"),
       subtitle = T, pValCutoff = 0.01, saveGeneList = T, addLabs = "KCNK16"
      )


### Fig 4f - gsea of the DGE results within Enterocyte_2 cluster
genes.df <- read.csv("./output/duod/linDEG/duod_Enterocyte_2geneList.csv")
geneListUp <- genes.df %>% arrange(p_val_adj) %>% filter(avg_log2FC > 0) %>% .$X
geneListDwn <- genes.df %>% arrange(p_val_adj) %>% filter(avg_log2FC < 0) %>% .$X

p <- plotGSEA(geneList = geneListUp, geneListDwn = geneListDwn, category = "C5", subcategory = NULL, size = 4.1,
              upCol = "blue", dwnCol = "red") + scale_x_continuous(limits = c(-4.5,6), name = "Signed log10(padj)") + theme( axis.title=element_text(size = 16),
                                                                                                                           title = element_text(size = 20),
                                                                                                                           plot.title = element_text(size = 20, hjust = 0.5)) + ggtitle("Enterocyte_2 gene ontology (CIE vs healthy)")
ggsave(paste("./output/", outName, "/", outName, "_fig4f.png", sep = ""), width = 9.5, height = 7)
