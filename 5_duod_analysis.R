#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")

#read in processed "All cells" dataset
seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/s3/230816_duod_h3c4_NoIntrons_res1.3_dims40_dist0.3_neigh50_S3.rds")
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
                  majorID ==  "epithelial" | majorID ==  "fibroblast")

table(seu.obj.sub$majorID)
table(seu.obj.sub$clusterID_2_1)
table(seu.obj.sub$orig.ident)


#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "./output/s2/", subName = "230916_duod_duod_h3c4_NoIntrons_2500", 
                     preSub = T, nfeatures = 2500, vars.to.regress = "percent.mt",  saveRDS = F
                    )

#clustree to determine clus resolution
# seu.obj <- readRDS(file = "./output/s2/230829_myeloid_duod_h3c4_NoIntrons_2000")
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


# seu.obj <- readRDS(file = "./output/s2/230717_myeloid_duod_h3c3_NoIntrons_2500_S2.rds")
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


ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "./output/cb_input/", 
               markers = paste0("./output/viln/",outName,"/",saveName,"_gene_list.csv"), 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", "majorID",
                                                   "clusterID_sub", "name2", "majorID_sub", "cellSource", "clusterID_final"), 
               skipEXPR = T,test = F,
               feats = c("SLC15A1", "ACE2", "IDO1", "CFTR", 
                         "BEST4", "NOTCH2", "AREG", "FABP3", 
                         "IRAG2", "SYNE2", "ONECUT2","SI")
                          
                          )



### Fig : Create UMAP by clusterID_sub
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              pt.size = 0.25,
            #  cols = majorColors.df$colz,
              label = T,
              label.box = T,
              shuffle = TRUE
)
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8) + NoLegend() #, labCol = majorColors.df$labCol
ggsave(paste("./output/", outName,"/", outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)


### Fig : Create UMAP by majorID_sub
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
ggsave(paste("./output/", outName,"/",outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)


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
ggsave(paste0("./output/", outName, "/", outName, "_rawUMAP_tiny.png"), width = 7, height = 7)



features <- c("SLC15A1", "SI", "ACE2", "TMEM37",
              #"MGAM","SLC7A8","PTPRH","CNDP1","APOC3","MICAL3","MGST3","TPT1","COX4I1",
              "ISG15","HERC6",
              "IDO1", "CFTR", "PLCB4","BEST4","NOTCH2","ADGRF5","CEACAM20",
              "DCN","COL1A1","CALD1",
              "AREG", "SYTL2", "AGR2", "FABP3", "FOXA3", "NRTN",
              "PTGS1", "IRAG2", "FYB1", "ALOX5", "TSNARE1", "SYNE2",
              "SCG2", "SCG3", "SCG5", "RFX6", "PCSK1", "RIMBP2")


#plot features
p <- majorDot(seu.obj = seu.obj, groupBy = "majorID_sub",
              yAxis = levels(seu.obj$majorID_sub),
                  features = features
                 ) + theme(axis.title = element_blank(),
                           axis.text = element_text(size = 12))
ggsave(paste("./output/", outName, "/", outName, "_majorDot.png", sep = ""), width =9, height = 4)



### Fig: Evlauate cell frequency by cluster
freqy <- freqPlots(seu.obj, method = 1, nrow= 2, groupBy = "majorID_sub", legTitle = "Cell source",refVal = "name2", showPval = T,
                   namez = unique(seu.obj$name2), 
                   colz = unique(seu.obj$colz)
                  ) + theme(strip.text = element_text(size=8))

res.ftest <- lapply(levels(freqy$data$majorID_sub), function (x){
    data.df <- freqy$data[freqy$data$majorID_sub ==  x, ]
    car::leveneTest(freq ~ cellSource, data = data.df)
})
              

ggsave(paste("./output/", outName, "/",outName, "_freqPlots.png", sep = ""), width = 8.5, height = 5)






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
ggsave(paste("./output/", outName, "/", outName, "_cieVSnorm_volcPlot.png", sep = ""), width = 7, height = 7)




###
createPB(seu.obj = seu.obj, groupBy = "majorID_sub", comp = "cellSource", biologicalRep = "name2", lowFilter = T, dwnSam =F, min.cell = 15,
                     clusters = "Enterocyte", outDir = paste0("./output/", outName,"/pseudoBulk/") , grepTerm = "H", grepLabel = c("Healthy","CIE") #improve - fix this so it is more functional
                    )

p <- pseudoDEG(metaPWD = paste0("./output/", outName,"/pseudoBulk/majorID_sub_deg_metaData.csv"), returnDDS = F, 
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("./output/", outName,"/pseudoBulk/"), outName = "allCells", idents.1_NAME = "CIE", idents.2_NAME = "Healthy",
          inDir = paste0("./output/", outName,"/pseudoBulk/"), title = "All cells", fromFile = T, meta = NULL, pbj = NULL, returnVolc = T, paired = F, pairBy = "", 
          minimalOuts = F, saveSigRes = T, filterTerm = "^ENSCAF", addLabs = NULL, mkDir = T
                     )

p  <- prettyVolc(plot = p[[2]], rightLab = NULL, leftLab = NULL, arrowz = F) + labs(x = "log2(FC) CIE vs Normal")
ggsave(paste("./output/", outName, "/", outName, "_cieVSnorm_volcPlot.png", sep = ""), width = 7, height = 7)


seu.obj$allCells <- "DGE analysis of epithelial cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
linDEG(seu.obj = seu.obj, groupBy = "allCells", comparision = "cellSource", outDir = paste0("./output/", outName,"/"), 
       outName = "epithelial", labCutoff = 10, contrast = c("CIE", "Healthy"),
       subtitle = T, pValCutoff = 0.01, saveGeneList = T, addLabs = "KCNK16"
      )

df.res <- read.csv("./output/duod/epithelial_DGE analysis of epithelial cellsgeneList.csv")
df.res %>% filter(avg_log2FC > 0, p_val < 0.01) %>% dim()
df.res %>% filter(avg_log2FC < 0, p_val < 0.01) %>% dim()



linDEG(seu.obj = seu.obj, groupBy = "majorID_sub", comparision = "cellSource", outDir = paste0("./output/", outName,"/"), 
       outName = "epithelial", labCutoff = 10, contrast = c("CIE", "Healthy"),
       subtitle = T, pValCutoff = 0.01, saveGeneList = T, addLabs = "KCNK16"
      )





files <- list.files(path = "./output/duod/linDEG/", pattern=".csv", all.files=FALSE,
                        full.names=T)

df.list <- lapply(files, read.csv, header = T)
df <- bind_rows(df.list, .id = "column_label")
df.res <- df[df$p_val_adj < 0.01, ] %>% group_by(cellType) %>% summarize(cnts = n()) %>% as.data.frame()
df.res <- df.res[match(levels(seu.obj$majorID_sub), df.res$cellType), ]
df.res$cnts <- as.numeric(df.res$cnts)

df.res <- df.res[!is.na(df.res$cellType), ]

rownames(df.res) <- df.res$cellType
df.res$cellType <- NULL

mat.res <- as.matrix(df.res)

png(file = paste0("./output/", outName, "/", outName, "deg_heatmap.png"), width=1250, height=2000, res=400)
par(mfcol=c(1,1))         

ht <- Heatmap(as.matrix(mat.res),#name = "mat", #col = col_fun,
              name = "# of DEGs",
              cluster_rows = F,
              row_title = "Cell type",
              col=viridis(100),
              cluster_columns = F,
              column_title = "# of DEGs",
              show_column_names = FALSE,
              column_title_side = "top",
              column_names_rot = 45,
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                                          labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm")),
              
              cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.0f", as.matrix(mat.res)[i, j]), x, y, gp = gpar(fontsize = 10, col = c(rep("black",7),rep("white",4))))
              })

draw(ht, padding = unit(c(2, 12, 2, 5), "mm"),show_heatmap_legend = FALSE)

dev.off()




lapply(levels(seu.obj$majorID_sub)[c(1,2,3)], function(x){
    dgea.df <- read.csv(paste0("./output/duod/linDEG/duod_",x,"geneList.csv"))
    geneListUp <- dgea.df %>% arrange(p_val_adj) %>% filter(avg_log2FC > 0) %>% .$X
    geneListDwn <- dgea.df %>% arrange(p_val_adj) %>% filter(avg_log2FC < 0) %>% .$X


    p <- plotGSEA(geneList = geneListUp, geneListDwn = geneListDwn, category = "C2")
    ggsave(paste("./output/", outName, "/linDEG/", outName, "_",x,"_enriched_terms_1.png", sep = ""), width = 9, height =7)

    p <- plotGSEA(geneList = geneListUp, geneListDwn = geneListDwn, category = "C5")
    ggsave(paste("./output/", outName, "/linDEG/", outName, "_",x,"_enriched_terms_2.png", sep = ""), width = 9, height =7)

    p <- plotGSEA(geneList = geneListUp, geneListDwn = geneListDwn, category = "C2",subcategory = "CP:REACTOME")
    ggsave(paste("./output/", outName, "/linDEG/", outName, "_",x,"_enriched_terms_3.png", sep = ""), width = 9, height =7)
})


### Fig : gsea of the DGE results
genes.df <- read.csv("./output/duod/linDEG/duod_Enterocyte_2geneList.csv")
geneListUp <- genes.df %>% arrange(p_val_adj) %>% filter(avg_log2FC > 0) %>% .$X
geneListDwn <- genes.df %>% arrange(p_val_adj) %>% filter(avg_log2FC < 0) %>% .$X

p <- plotGSEA(geneList = geneListUp, geneListDwn = geneListDwn, category = "C5", subcategory = NULL,
              upCol = "blue", dwnCol = "red") + scale_x_continuous(limits = c(-4.5,6), name = "Signed log10(padj)") + theme( axis.title=element_text(size = 16),
                                                                                                                           title = element_text(size = 20),
                                                                                                                           plot.title = element_text(size = 20, hjust = 0.5)) + ggtitle("Enterocyte_2 gene ontology (CIE vs healthy)")
ggsave(paste("./output/", outName, "/", outName, "_enriched_terms_duod_Enterocyte_2.png", sep = ""), width = 9.5, height = 7)



viln.df <- read.csv("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/viln/duod/231022_duod_duod_h3c4_NoIntrons_2500_gene_list.csv")
viln.df[viln.df$cluster == "Progenitor", ] %>% top_n(-24, wt = p_val_adj) %>% pull(gene) %>% paste0(., collapse = ",")


seu.gut.duod <- MuDataSeurat::ReadH5AD("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/epi_log_counts02_v2.h5ad")
reference <- seu.gut.duod

seu.gut.Mesenchyme <- MuDataSeurat::ReadH5AD("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/Mesenchyme_log_counts02_v2.h5ad")
reference <- seu.gut.Mesenchyme

reference <- SCTransform(reference, verbose = T)
reference <- RunPCA(reference)
# DefaultAssay(reference) <- "integrated"

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

pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "predicted.id",
              #cols = levels(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.25,
              split.by = "predicted.id",
              ncol = 4,
              label = F,
              label.box = F,
              shuffle = F
)
pi <- formatUMAP(plot = pi)
ggsave(paste("./output/", outName, "/",outName, "_umap_Predicted.png", sep = ""), width = 14, height = 14)




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
ggsave(paste("./output/", outName, "/",outName, "splitFeats.png", sep = ""), width = 2, height = 4)
