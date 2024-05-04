#!/usr/bin/Rscript

#load custom functions & packages
# source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")

### Analysis note: 
# This script loads in the previously processed Seurat object (./output/s3/230816_duod_h3c4_NoIntrons_res1.3_dims40_dist0.3_neigh50_S3.rds)
# then subsets on epithelial cells and generates all figures assocaited with Figure 4

################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin epithelial preprocessing   ######## <<<<<<<<<<<<<<
################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#read in processed "All cells" dataset
seu.obj <- readRDS("./output/s3/230816_duod_h3c4_NoIntrons_res1.3_dims40_dist0.3_neigh50_S3.rds")
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
                  majorID ==  "epithelial")

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

### Note cluster 10 looks to be T cells -- remove and re-integrate
seu.obj.sub <- subset(seu.obj, invert = T,
                      subset = 
                      clusterID_sub ==  "10")
table(seu.obj.sub$clusterID_sub)
table(seu.obj.sub$orig.ident)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "./output/s2/", subName = "230916_duod_duod_h3c4_NoIntrons_2500", 
                     preSub = T, nfeatures = 2500, vars.to.regress = "percent.mt",  saveRDS = F
                    )

#run clustree
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "230916_duod_duod_h3c4_NoIntrons_2500", test_dims = c("40","35", "30"), algorithm = 3, prefix = "integrated_snn_res.")

#visulize the data and proceed with analysis
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230916_duod_duod_h3c4_NoIntrons_2500", final.dims = 40,
                       final.res = 1.7, stashID = "clusterID_sub", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.25,
                       n.neighbors = 25, assay = "integrated", saveRDS = T,
                       features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


############################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin epithelial analysis   ######## <<<<<<<<<<<<<<
############################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in preprocessed data
seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/s3/230916_duod_duod_h3c4_NoIntrons_2500_res1.7_dims40_dist0.25_neigh25_S3.rds")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))
colz.df <- read.csv("./cellColz.csv", header = F)
colz.df <- colz.df[colz.df$V2 == "duod", ]
outName <- "duod"

#stash the new idents -- note: had to overcluster the data to detect tuft, stromal, and enteroenddocrine, so this collapses the enterocytes base on the clustree output
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


### supp data - cell type gene signatures

#generate by clusterID (note: this was overclustered)
# vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = "230916_duod_duod_h3c4_NoIntrons_2500",
#                      outDir = "./output/viln/duod/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
#                     )

#generate based on cell types
saveName <- "240128_duod_duod_h3c4_NoIntrons_2500"
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_sub", numOfFeats = 24, outName = saveName,
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )

#append with surface call
surface.markers <- read.csv("./surface_master.csv")[ ,c("UniProt.gene", "UniProt.description", "Surfaceome.Label", "Surfaceome.Label.Source")] %>% filter(!duplicated(UniProt.gene))
cluster.markers <- read.csv("./output/viln/duod/240128_duod_duod_h3c4_NoIntrons_2500_gene_list.csv")
write.csv(cluster.markers[ ,c(7,8,2:6)] %>% left_join(surface.markers, by = c("gene" = "UniProt.gene")),
          file = "./output/supplementalData/supplemental_data_8.csv", row.names = F)


### supp data - export data fpr cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "./output/cb_input/", 
               markers = paste0("./output/viln/",outName,"/",saveName,"_gene_list.csv"), 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", "majorID",
                                                   "clusterID_sub", "name2", "majorID_sub", "cellSource", "clusterID_final"), 
               skipEXPR = T,test = F,
               feats = c("SLC15A1", "ACE2", "IDO1", "CFTR", 
                         "BEST4", "NOTCH2", "AREG", "FABP3", 
                         "IRAG2", "SYNE2", "ONECUT2","SI")
                          
                          )


### Fig supp 6a - UMAP by clusterID_sub
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = TRUE
) + NoLegend()
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, smallAxes = T) #, labCol = majorColors.df$labCol
ggsave(paste("./output/", outName,"/", outName, "_supp6a.png", sep = ""), width = 7, height = 7)


### Fig 4a - UMAP by majorID_sub
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_final",
              pt.size = 0.25,
              cols = colz.df$V1,
              label = F,
              label.box = F,
              shuffle = TRUE
) + NoLegend() 
p <- formatUMAP(plot = pi, smallAxes = T)
ggsave(paste0("./output/", outName, "/", outName, "_fig4a.png"), width = 7, height = 7)


### Fig 4a: create labels to be cropped onto UMAP
majorColors.df <- as.data.frame(rep("", length(colz.df$V1)))
colnames(majorColors.df) <- "clusterID_final"
majorColors.df$colz <- colz.df$V1
majorColors.df$title <- ""
majorColors.df$labCol <- "black"
majorColors.df$labz <- ( 1:length(colz.df$V1) - 1)
leg <- cusLeg(legend = majorColors.df, clusLabel = "labz", legLabel = "clusterID_final", colorz = "colz",labCol = "labCol",colz = 1, rowz = NULL, groupLabel = "title", dotSize = 10, groupBy = "title",sortBy = "labz", compress_x = 0.9, textSize = 6)
ggsave(paste("./output/", outName, "/", outName, "_leg_forUMAP_labels.png", sep = ""), width = 2.75, height = 7)


### Fig 4b - dotplot of major featrures
features <- c("SLC15A1", "SI", "ACE2", "TMEM37",
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


features <- c("SI","SLC22A4","SLC43A2", "SLC40A1",
              "SLC39A4", "SLC25A5", "SLC25A6", "SLC51B",
              "ACTA2", "BEST4", "CFTR", "IFIT1")
p2 <- prettyFeats(seu.obj = seu.obj, nrow = 3, ncol = 4, features = features, 
             color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", outName, "_fig4b.png", sep = ""), width = 12, height = 9)

### Fig supp 6b - reference map using the human gut atlas -- epithieal reference
seu.gut.duod <- MuDataSeurat::ReadH5AD("./epi_log_counts02_v2.h5ad")
#download reference from https://www.gutcellatlas.org/#datasets

#prep the ref
reference <- seu.gut.duod
reference <- SCTransform(reference, verbose = FALSE)
reference <- RunPCA(reference)

#map the data
anchors <- FindTransferAnchors(
    reference = reference,
    query = seu.obj,
    normalization.method = "SCT",
    reference.reduction = "pca", #reference.reduction = "umap",
    dims= 1:50 #dims= 1:2
)

#transfer annotations
predictions <- TransferData(anchorset = anchors, refdata = reference$annotation,
    dims = 1:50)
seu.obj <- AddMetaData(seu.obj, metadata = predictions)
seu.obj$predicted.human <- seu.obj$predicted.id

#plot the data
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
ggsave(paste("./output/", outName,"/",outName, "_sup6b.png", sep = ""), width = 7, height = 7)


### Fig supp 6c - reference map using the human gut atlas -- mesnechymal reference
reference <- MuDataSeurat::ReadH5AD("./Mesenchyme_log_counts02_v2.h5ad")
#download reference from https://www.gutcellatlas.org/#datasets

#prep the ref
reference <- SCTransform(reference, verbose = FALSE)
reference <- RunPCA(reference)

#complete mapping
anchors <- FindTransferAnchors(
    reference = reference,
    query = seu.obj,
    normalization.method = "SCT",
    reference.reduction = "pca", #reference.reduction = "umap",
    dims= 1:50 #dims= 1:2
)

#transfer labels
predictions <- TransferData(anchorset = anchors, refdata = reference$annotation,
    dims = 1:50)
seu.obj <- AddMetaData(seu.obj, metadata = predictions)
seu.obj$predicted.human.Mesenchyme <- seu.obj$predicted.id

#plot the data
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
ggsave(paste("./output/", outName,"/",outName, "_sup6c.png", sep = ""), width = 7, height = 7)


### Fig sup 6d: plot enrichment scores
ecLists <- read.csv("./gut_ecTerms.csv", header = T)
#download reference from https://www.gutcellatlas.org/#datasets
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
ggsave(paste("./output/", outName, "/", outName, "_sup6d.png", sep = ""), width = 8, height = 6)


### Fig 4c - evaluate cell frequency by cluster
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

#get data for in text call
freqy$data %>% filter(majorID_sub == "Tuft") %>% group_by(cellSource) %>% summarize(ave = mean(freq))

log2FD_threshold <- 1

Idents(seu.obj) <- "name2"
set.seed(12)
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$name2)))
prop_test <- sc_utils(seu.obj.sub)
prop_test <- permutation_test( prop_test, cluster_identity = "majorID_sub", sample_1 = "Healthy", sample_2 = "CIE", sample_identity = "cellSource" )
p <- permutation_plot(prop_test)  + theme(axis.title.y = element_blank(),
                                          legend.position = "top") + guides(colour = guide_legend("", nrow = 2, 
                                                                                                  byrow = TRUE)) + coord_flip()
res.df <- prop_test@results$permutation
res.df <- res.df %>% mutate(Significance = as.factor(ifelse(obs_log2FD < -log2FD_threshold & FDR < 0.01 & boot_CI_97.5 < -log2FD_threshold,"Down",
                                                            ifelse(obs_log2FD > log2FD_threshold & FDR < 0.01 & boot_CI_2.5 > log2FD_threshold, "Up", "n.s.")))
                           ) %>% arrange(obs_log2FD)
res.df$clusters <- factor(res.df$clusters, levels = c(res.df$clusters))

p <- ggplot(res.df, aes(x = clusters, y = obs_log2FD)) + 
geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, 
                    color = Significance)) + theme_bw() + geom_hline(yintercept = log2FD_threshold, 
                                                                     lty = 2) + geom_hline(yintercept = -log2FD_threshold, 
                                                                                           lty = 2) + 
geom_hline(yintercept = 0) + scale_color_manual(values = c("blue", "grey")
                                               ) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
                                                         axis.title.x = element_blank(),
                                                         legend.position = "top",
                                                         plot.margin = margin(c(3,3,0,10))
                                                        ) + ylab("abundance change (log2FC)")
ggsave(paste("./output/", outName, "/",outName, "_fig2c.png", sep = ""), width = 3.5, height = 2, scale = 2 )


### Fig supp 6e - pseudobulk dge analysis
seu.obj$allCells <- "Epithelial cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
createPB(seu.obj = seu.obj, groupBy = "allCells", comp = "cellSource", biologicalRep = "name2", lowFilter = T, dwnSam =F, min.cell = 15,
                     clusters = NULL, outDir = paste0("./output/", outName,"/pseudoBulk/") , grepTerm = "H", grepLabel = c("Healthy","CIE") #improve - fix this so it is more functional
                    )
p <- pseudoDEG(metaPWD = paste0("./output/", outName,"/pseudoBulk/allCells_deg_metaData.csv"), returnDDS = F, 
          padj_cutoff = 0.05, lfcCut = 1, outDir = paste0("./output/", outName,"/pseudoBulk/"), outName = "allCells", idents.1_NAME = "CIE", idents.2_NAME = "Healthy",
          inDir = paste0("./output/", outName,"/pseudoBulk/"), title = "All cells", fromFile = T, meta = NULL, pbj = NULL, returnVolc = T, paired = F, pairBy = "", 
          minimalOuts = F, saveSigRes = T, filterTerm = "^ENSCAF", addLabs = NULL, mkDir = T, labSize = 5.5, strict_lfc = F,
                     )
pi  <- prettyVolc(plot = p[[1]], rightLab = NULL, leftLab = NULL, arrowz = F) + labs(title = "CIE vs Healthy (within all myeloid cells)", x = "log2(Fold change)") + NoLegend() + theme(panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 2),
                                      axis.line = element_blank(),
                                      plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
ggsave(paste("./output/", outName, "/", outName, "_supp6e.png", sep = ""), width = 7, height = 7)


### Fig 4d - dge scatter plot
seu.obj$allCells <- "DGE analysis of epithelial cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
linDEG(seu.obj = seu.obj, groupBy = "allCells", comparision = "cellSource", outDir = paste0("./output/", outName,"/fig4d_"), 
       outName = "epithelial", labCutoff = 10, contrast = c("CIE", "Healthy"),
       subtitle = T, pValCutoff = 0.01, logfc.threshold = 0.58, saveGeneList = T, addLabs = "KCNK16"
      )

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


### Fig extra - DEGs by samples
lapply(unique(seu.obj$name2), function(sampleName){
    seu.obj.sub <- subset(seu.obj, subset = name2 == sampleName)
    
    pi <- DimPlot(seu.obj.sub, 
                  reduction = "umap", 
                  group.by = "majorID_sub",
                  pt.size = 0.25,
                  label = F,
                  label.box = F,
                  repel = F,
                  cols = colz.df$V1
    )
    p1 <- formatUMAP(plot = pi) + theme(plot.title = element_text(size= 20)) + ggtitle(unique(seu.obj.sub$name2))
    
#     features <- c("IL1B", "S100A12", 
#                   "SOD2", "KCNK16")
    
    features <- c("KCNK16", "ANKRD37", 
                  "CYP1B1", "HSPA2")
    p2 <- prettyFeats(seu.obj = seu.obj.sub, nrow = 2, ncol = 2, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)

    p <- p1 + p2
    ggsave(plot = p, paste0("./output/", outName, "/", unique(seu.obj.sub$name2), "_key_feats.png"), width = 12, height = 6, device = "png")
})


### Fig extra - dge scatter plot within each epithelial cell type
linDEG(seu.obj = seu.obj, groupBy = "majorID_sub", comparision = "cellSource", outDir = paste0("./output/", outName,"/"), 
       outName = "epithelial", labCutoff = 10, contrast = c("CIE", "Healthy"),
       subtitle = T, pValCutoff = 0.01, logfc.threshold = 0.58, saveGeneList = T, addLabs = "KCNK16"
      )


### Fig supp 6e - pseudobulk dge analysis
createPB(seu.obj = seu.obj, groupBy = "majorID_sub", comp = "cellSource", biologicalRep = "name2", lowFilter = T, dwnSam =F, min.cell = 15,
                     clusters = NULL, outDir = paste0("./output/", outName,"/pseudoBulk/") , grepTerm = "H", grepLabel = c("Healthy","CIE") #improve - fix this so it is more functional
                    )
p <- pseudoDEG(metaPWD = paste0("./output/", outName,"/pseudoBulk/majorID_sub_deg_metaData.csv"), returnDDS = F, 
          padj_cutoff = 0.1, lfcCut = 1, outDir = paste0("./output/", outName,"/pseudoBulk/"), outName = "allCells", idents.1_NAME = "CIE", idents.2_NAME = "Healthy",
          inDir = paste0("./output/", outName,"/pseudoBulk/"), title = "All cells", fromFile = T, meta = NULL, pbj = NULL, returnVolc = T, paired = F, pairBy = "", 
          minimalOuts = F, saveSigRes = T, filterTerm = "^ENSCAF", addLabs = NULL, mkDir = T, labSize = 5.5, strict_lfc = F,
                     )


### Supp fig xx -- heatmap of degs by each cluster
files <- lapply(levels(seu.obj$majorID_sub), function(x){
    list.files(path = paste0("./output/", outName,"/pseudoBulk/", x), pattern = ".csv", all.files = FALSE, full.names = T)
})

df.list <- lapply(unlist(files), read.csv, header = T)

seu.obj$type <- factor(paste0(as.character(seu.obj$majorID_sub), "--", as.character(seu.obj$cellSource)),
                       levels = paste0(rep(levels(seu.obj$majorID_sub), each = 2), "--", c("Healthy", "CIE")))

res.df <- do.call(rbind, df.list) %>% 
    filter(abs(log2FoldChange) > 1) %>% 
    filter(!grepl("^ENS", gene)) %>%
    group_by(gs_base) %>%
    top_n(10, wt = abs(log2FoldChange))

sig.mat <- matrix(nrow = length(unique(res.df$gene)), ncol = length(levels(seu.obj$type)),
                  dimnames = list(unique(res.df$gene),
                                  toupper(levels(seu.obj$type))))

for(i in 1:nrow(sig.mat)){
    for(j in 1:ncol(sig.mat)){
        cellType <- strsplit(colnames(sig.mat)[j], "--")[[1]][1]
        condition <- strsplit(colnames(sig.mat)[j], "--")[[1]][2]
        if(cellType %in% res.df[res.df$gene == rownames(sig.mat)[i], ]$gs_base){
            lfc <- res.df[res.df$gene == rownames(sig.mat)[i] & res.df$gs_base == cellType, ]$log2FoldChange
            if(lfc > 1 & condition == "CIE"){
                sig.mat[i, j] <- "*"
            } else if(lfc < -1 & condition == "HEALTHY"){
                sig.mat[i, j] <- "*"
            } else{
                sig.mat[i, j] <- ""
            }
        } else{
            sig.mat[i, j] <- ""
        }
    }
}

res.df <- res.df[!duplicated(res.df$gene), ]


#extract metadata and data
metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$RNA@counts)) #use raw count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type

#get cell type expression averages - do clus avg expression by sample
clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL

#filter matrix for DEGs and scale by row
clusAvg_expression <- clusAvg_expression[ ,colnames(clusAvg_expression) %in% res.df$gene]
mat_scaled <- t(apply(t(log1p(clusAvg_expression)), 1, scale))
colnames(mat_scaled) <- rownames(clusAvg_expression)
mat_scaled <- mat_scaled[ ,match(colnames(sig.mat), toupper(colnames(mat_scaled)))]
mat_scaled <- mat_scaled[match(rownames(sig.mat), rownames(mat_scaled)), ]  

#set annotations
clus <- colz.df$V1
names(clus) <- levels(seu.obj$majorID_sub)
cond_colz <- c("mediumseagreen","mediumpurple1")
names(cond_colz) <- c("Healthy","CIE")

# heat_col <- viridis(option = "magma",100)
ha <- HeatmapAnnotation(
    Cluster = factor(unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"--")[[1]][1]})),
                     levels = levels(seu.obj$majorID_sub)),
    Condition = unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"--")[[1]][2]})),
    border = TRUE,
    col = list(Cluster = clus, Condition = cond_colz),
    show_annotation_name = FALSE
)

#plot the data
ht <- Heatmap(
    mat_scaled,
    name = "mat",
    cluster_rows = F,
    row_names_gp = gpar(fontsize = 12),
    show_row_names = T,
    cluster_columns = F,
    top_annotation = ha,
    show_column_names = F,
    column_split = factor(unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"--")[[1]][1]})),
                          levels = levels(seu.obj$majorID_sub)),
    row_title = NULL,
    column_title = NULL,
    heatmap_legend_param = list(
        title = "Scaled expression",
        direction = "horizontal"
        ),
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sig.mat[i, j], x, y, gp = gpar(fontsize = 12, col = "black"))
    }
)

png(file = paste0("./output/", outName, "/", outName, "_fig3e.png"), width=2500, height=3000, res=400)
par(mfcol=c(1,1))   
draw(ht, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "bottom")

for(i in 1:length(levels(seu.obj$majorID_sub))){
    decorate_annotation("Cluster", slice = i, {
        grid.text(paste0("c", (1:length(levels(seu.obj$majorID_sub))) - 1)[i], just = "center")
    })
}
dev.off()


### Fig 4f - gsea of the DGE results within Enterocyte_2 cluster
genes.df <- do.call(rbind, df.list) %>% 
    filter(abs(log2FoldChange) > 1) %>% 
    filter(!grepl("^ENS", gene)) %>%
    filter(grepl("ENTEROCYTE", gs_base))
# genes.df <- read.csv("./output/duod/pseudoBulk/Enterocyte_2/allCells_cluster_Enterocyte_2_all_genes.csv")
geneListUp <- genes.df %>% arrange(padj) %>% filter(log2FoldChange > 1) %>% pull(gene)
geneListDwn <- genes.df %>% arrange(padj) %>% filter(log2FoldChange < -1) %>% pull(gene)

p <- plotGSEA(geneList = geneListUp, geneListDwn = geneListDwn, category = "C5",
              subcategory = NULL, size = 3.5, termsTOplot = 16, upCol = "blue", 
              dwnCol = "red", lolli = T) + 
    scale_x_continuous(
        limits = c(-10,10), 
        name = "Signed log10(padj)", 
        breaks = c(-ceiling(max(abs(p$data$x_axis))*1.05), -ceiling(max(abs(p$data$x_axis))*1.05)/2, 0, 
                   ceiling(max(abs(p$data$x_axis))*1.05)/2, ceiling(max(abs(p$data$x_axis))*1.05))
    ) + 
    theme(
        axis.title=element_text(size = 16),
        title = element_text(size = 20),
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.margin = margin(7, 7, 7, 42, "pt"),
        legend.position = c(0.9, 0.9),
    ) + 
    ggtitle("Enterocyte gene ontology\n(CIE vs healthy)")
ggsave(paste("./output/", outName, "/", outName, "_fig4f.png", sep = ""), width = 10, height = 8)

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end epithelial analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<