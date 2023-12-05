#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")

### complete analysis with n=4 cie

#################################################
##### Analysis of n3 healthy & n4 CIE douds #####
#################################################

### Metadata transfer from indivudally annotated datasets

#load in processed allCells data
seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/s3/230816_duod_h3c4_NoIntrons_res1.3_dims40_dist0.3_neigh50_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./colorID_cie3v4.csv", groupBy = "clusterID_2_1", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name2")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))
outName <- "allCells"
subname <- "n3n4"


Idents(seu.obj) <- "majorID"
seu.obj <- RenameIdents(seu.obj, c("tcell" = "T cell", "epithelial" = "Epithelial", 
                                   "myeloid" = "Myeloid", "bcell" = "B cell",
                                   "cycling" = "Cycling T cell", "mast" = "Mast cell",
                                   "plasma" = "Plasma cell")
                       )


seu.obj$majorID_pertyName <- Idents(seu.obj)

seu.obj <- subset(seu.obj,
                  subset = majorID_pertyName == "Cycling T cell" |  majorID_pertyName == "B cell" |  majorID_pertyName == "Mast cell" | majorID_pertyName == "Plasma cell")

table(seu.obj$majorID_pertyName)

ct.l3 <- droplevels(seu.obj$majorID_pertyName)

#load in processed epitheilal subset data
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

ct.l3 <- c(ct.l3,seu.obj$majorID_sub)

#load in processed myeloid subset data
seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/s3/230829_myeloid_duod_h3c4_NoIntrons_2500_res0.5_dims40_dist0.25_neigh25_S3.rds")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))
colz.df <- read.csv("./cellColz.csv", header = F)
colz.df <- colz.df[colz.df$V2 == "myeloid", ]
outName <- "myeloid"


Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "Monocyte", "1" = "Neutrophil", 
                                   "2" = "Eosinophil", "3" = "Macrophage",
                                   "4" = "Neutrophil", "5" = "Eosinophil",
                                   "6" = "cDC1", "7" = "IL22RA2_DC",
                                   "8" = "Eosinophil")
                       )

seu.obj$majorID_sub <- Idents(seu.obj)
ct.l3 <- c(ct.l3,seu.obj$majorID_sub)

#load in processed tcell subset data
seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/s3/230913_tcell_duod_h3c4_NoIntrons_2500_res0.6_dims35_dist0.3_neigh30_S3.rds")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))
colz.df <- read.csv("./cellColz.csv", header = F)
outName <- "tcell"

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "Trm", "1" = "Trm", 
                                   "2" = "Tnaive", "3" = "Trm",
                                   "4" = "gdT_1", "5" = "Tnaive",
                                   "6" = "NK_T", "7" = "Tnaive",
                                   "8" = "Tnaive", "9" = "gdT_2",
                                   "10" = "Trm", "11" = "Treg",
                                   "12" = "Trm", "13" = "ILC2",
                                   "14" = "NK","15" = "T_IFN")
                       )


seu.obj$majorID_sub <- Idents(seu.obj)
ct.l3 <- c(ct.l3,seu.obj$majorID_sub)

#load in processed allCells data again... to add new metadata
seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/s3/230816_duod_h3c4_NoIntrons_res1.3_dims40_dist0.3_neigh50_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./colorID_cie3v4.csv", groupBy = "clusterID_2_1", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name2")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))
outName <- "allCells"
subname <- "n3n4"

Idents(seu.obj) <- "majorID"
seu.obj <- RenameIdents(seu.obj, c("tcell" = "T cell", "epithelial" = "Epithelial", 
                                   "myeloid" = "Myeloid", "bcell" = "B cell",
                                   "cycling" = "Cycling T cell", "mast" = "Mast cell", "plasma" = "Plasma cell")
                       )


seu.obj$majorID_pertyName <- Idents(seu.obj)


seu.obj <- AddMetaData(seu.obj, metadata = ct.l3, col.name = "celltype.l3")

seu.obj$celltype.l3 <- as.factor(ifelse(is.na(seu.obj$celltype.l3), "remove", as.character(seu.obj$celltype.l3)))
seu.obj <- subset(seu.obj, invert = T,
                  subset = celltype.l3 == "remove"
                 )

seu.obj$celltype.l3 <- droplevels(seu.obj$celltype.l3)
# saveRDS(seu.obj, "./output/s3/canine_duodenum_annotated.rds")


### Proper analysis of the dataset

#load in annotated dataset
seu.obj <- readRDS("./output/s3/canine_duodenum_annotated.rds")
outName <- "allCells"
subname <- "n3n4"

clusterID_final <- table(seu.obj$majorID_pertyName) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
mutate(clusterID_final=row_number()-1) %>% arrange(clusterID_final) 

newID <- clusterID_final$clusterID_final
names(newID) <- clusterID_final$Var1
Idents(seu.obj) <- "majorID_pertyName"
seu.obj <- RenameIdents(seu.obj, newID)
table(Idents(seu.obj))
seu.obj$clusterID_final <- Idents(seu.obj)


### Fig extra - check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


### Supp data - run FindAllMarkers on major and high-res cell types
saveName <- "231022_allCells_duod_h3c4_NoIntrons_2500"
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_pertyName", numOfFeats = 24, outName = saveName,
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )

saveName <- "231204_allCells_duod_h3c4_NoIntrons_2500"
vilnPlots(seu.obj = seu.obj, groupBy = "celltype.l3", numOfFeats = 24, outName = saveName,
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )


### Supp data - export data for cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "./output/cb_input/", 
               markers = paste0("./output/viln/",outName,"/",saveName,"_gene_list.csv"), 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", "majorID",
                                                   "clusterID_sub", "name2", "majorID_pertyName", "celltype.l3", "cellSource", "clusterID_final"), 
               skipEXPR = F,test = F,
               feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                         "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                         "CD4", "MS4A1", "RELN","EPCAM")
                          
                          )


### Fig supp 1a - create UMAP by celltype.l3
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "celltype.l3",
              pt.size = 0.25,
            #  cols = colz.df$V1,
              label = T,
              label.box = T,
              shuffle = TRUE,
              repel = T
)
p <- formatUMAP(plot = pi) + NoLegend() + theme(plot.title = element_text(size = 18, vjust = 1),
                                                axis.title = element_blank(),
                                                panel.border = element_blank()) + ggtitle("Annotated canine duodenum atlas")

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
ggsave(paste("./output/", outName, "/", outName, "_ctl3_UMAP.png", sep = ""), width = 7, height = 7)


### Fig 1a - create UMAP by major cell types
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_final",
              cols = c("#FF89B3", "#C89504", "#00ABFD", "#B983FF", "#FA7476", "#0A9B9F", "#9B8EFF"),
              pt.size = 0.25,
              label = T,
              label.box = T
 )
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.8, labCol = "black") + NoLegend() + theme(axis.title = element_blank(),
                                               panel.border = element_blank(),
                                              plot.margin = unit(c(-7, -7, -7, -7), "pt"))

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

pi <- p + inset_element(axes,left= 0,
  bottom = 0,
  right = 0.25,
  top = 0.25,
                       align_to = "full")
ggsave(paste("./output/", outName,  "/", subname,"/", outName, "_fig1a.png", sep = ""), width = 7, height = 7)


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
p <- formatUMAP(plot = pi) + NoLegend() + theme(plot.title = element_text(size = 18, vjust = 1),
                                               axis.title = element_blank(),
                                               panel.border = element_blank(),
                                                plot.margin = unit(c(-7, -7, -7, -7), "pt")) + ggtitle("Canine PBMC reference mapping")
ggsave(paste("./output/", outName,"/",outName, "_umap_Predicted_PBMC.png", sep = ""), width = 7, height = 7)

gc()

### Fig - reference map using human epithelial data
reference <- MuDataSeurat::ReadH5AD("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/epi_log_counts02_v2.h5ad")

reference <- SCTransform(reference, verbose = FALSE)
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
seu.obj$predicted.human_epi <- seu.obj$predicted.id

pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "predicted.id",
              #cols = levels(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F
)
p <- formatUMAP(plot = pi) + NoLegend() + theme(plot.title = element_text(size = 18, vjust = 1),
                                               axis.title = element_blank(),
                                               panel.border = element_blank(),
                                                plot.margin = unit(c(-7, -7, -7, -7), "pt")) + ggtitle("Human epithelial reference mapping")
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_umap_Predicted_hu_epi.png", sep = ""), width = 7, height = 7)

gc()

### Fig - reference map using human epithelial data
reference <- MuDataSeurat::ReadH5AD("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/Tcell_log_counts02_v2.h5ad")

reference <- SCTransform(reference, verbose = FALSE)
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
seu.obj$predicted.human_mesen <- seu.obj$predicted.id

pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "predicted.id",
              #cols = levels(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F
)
p <- formatUMAP(plot = pi) + NoLegend() + theme(plot.title = element_text(size = 18, vjust = 1),
                                               axis.title = element_blank(),
                                               panel.border = element_blank(),
                                                plot.margin = unit(c(-7, -7, -7, -7), "pt")) + ggtitle("Human T cell reference mapping")
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_umap_Predicted_hu_tcell.png", sep = ""), width = 7, height = 7)

gc()


### Fig 1b - dot plot by major cell types
fig1c <- majorDot(seu.obj = seu.obj, groupBy = "majorID_pertyName",
                  yAxis = c("Epithelial","T cell","Myeloid","Plasma cell","B cell","Mast cell","Cycling T cell"),
                  features = c("SI", "FABP1", "RBP2", "GUCA2A", "APOA1", 
                               "CD3E", "CCL4", "GZMA","GZMB", "CD8A",
                               "AIF1", "C1QC", "S100A12", "LYZ", "CXCL8", 
                               "JCHAIN", "RARRES2", "IGHM", "MS4A1",
                                "PAX5","CD22","CD19","BANK1",
                               "KIT", "IGF1", "IL4I1",
                               "TOP2A","MKI67")
                 ) + theme(axis.title = element_blank(),
                           axis.text = element_text(size = 12)) +  scale_colour_viridis(option="magma", name='Average\nexpression', 
                                                                                        breaks = c(-0.5, 1, 2),
                                                                                        labels = c("-0.5", "1", "2")
                                                                                       ) + guides(color = guide_colorbar(title = 'Scaled\nExpression  ')) 
ggsave(paste("./output/", outName, "/", subname, "/", outName, "_fig1b.png", sep = ""), width = 8, height = 4)


### Fig 1c - feature plots of key features
features <- c("PTPRC", "CD3G", "IL7R", 
              "SI", "MUC13", "ANPEP", 
              "COL1A1", "AIF1", "LYZ",
              "JCHAIN", "TOP2A", "CENPF",
              "KIT", "PAX5", "MS4A1"
             )

titlez <- c("PTPRC (CD45)", "CD3G", "IL7R", 
              "SI", "MUC13", "ANPEP", 
              "COL1A1", "AIF1 (Iba1)", "LYZ",
              "JCHAIN", "TOP2A", "CENPF",
              "KIT", "PAX5", "MS4A1 (CD20)"
           )

colz <- c("black", "#FF89B3", "#FF89B3", 
          "#C89504", "#C89504", "#C89504",
          "#C89504", "#00ABFD", "#00ABFD",
          "#B983FF", "#FA7476", "#FA7476", 
          "#0A9B9F", "#9B8EFF", "#9B8EFF"
         )

p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, features = features, color = colz, order = T, titles = titlez, returnPlots = F, title.size = 16, pt.size = 0.00000001, min.cutoff = "q10", noLegend = T) 
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_fig1c.png", sep = ""), width = 9, height = 15)


### Fig 1d - stacked bar graph by major cell types
Idents(seu.obj) <- "cellSource"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data$cellSource)))
seu.obj.ds$majorID_pertyName <- factor(seu.obj.ds$majorID_pertyName, levels = rev(c("Epithelial","T cell","Myeloid","Plasma cell","B cell","Mast cell","Cycling T cell")))

p <- stackedBar(seu.obj = seu.obj.ds, downSampleBy = "cellSource", groupBy = "name2", clusters = "majorID_pertyName") +
scale_fill_manual(labels = levels(seu.obj$name2), 
                  values = levels(seu.obj$colz)) + theme(axis.title.y = element_blank(),
                                                         axis.title.x = element_text(size = 14),
                                                         axis.text = element_text(size = 12)
                                                        ) 
ggsave(paste("./output/", outName,"/", subname, "/",outName, "_fig1d.png", sep = ""), width =7, height = 4)


### Fig supp 1e: umap by sample
Idents(seu.obj) <- "cellSource"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data$cellSource)))
pi <- DimPlot(seu.obj.ds, 
              reduction = "umap", 
              group.by = "name2",
              cols = levels(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.25,
              label = FALSE,
              shuffle = TRUE
)
p <- formatUMAP(pi) + labs(colour="") + theme(axis.title = element_blank(),
                                               panel.border = element_blank(),
                                                plot.margin = unit(c(-7, -7, -7, -7), "pt"),
                                                  legend.position = "top", 
                                              legend.justification = "center",
                                              legend.direction = "horizontal",legend.title=element_text(size=12)) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))
ggsave(paste("./output/", outName, "/", subname, "/", outName, "_supp1e.png", sep = ""), width =7, height = 7)


### Fig supp 1e: stats by cie vs healthy
freqy <- freqPlots(seu.obj, method = 2, nrow= 2, groupBy = "majorID_pertyName", legTitle = "Cell source",refVal = "name2", showPval = T,
              namez = "name2", 
              colz = "colz"
              )

ggsave(paste("./output/", outName, "/", subname, "/",outName, "_supp1f.png", sep = ""), width = 8, height = 6)

#run levene test to ensure proper test was used
res.ftest <- lapply(levels(freqy$data$majorID_pertyName), function(x){
    data.df <- freqy$data[freqy$data$majorID_pertyName ==  x, ]
    leveneTest(freq ~ cellSource, data = data.df)
})
               

### Fig 1e: dge analysis all cells
seu.obj$allCells <- "DGE analysis of all cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
linDEG(seu.obj = seu.obj, threshold = 1, thresLine = F, groupBy = "allCells", comparision = "cellSource", contrast = c("CIE", "Healthy"),
       outDir = paste0("./output/", outName,"/", subname,"/"), 
       outName = "all_cells", cluster = NULL, labCutoff = 15, noTitle = F,
                   colUp = "red", colDwn = "blue", subtitle = T, returnUpList = F, returnDwnList = F, forceReturn = F, useLineThreshold = F, pValCutoff = 0.01, saveGeneList = T, addLabs = ""
                  )


### Fig extra: dge analysis by major cell types
linDEG(seu.obj = seu.obj, groupBy = "majorID_pertyName", comparision = "cellSource", contrast= c("CIE","Healthy"),
       outDir = paste0("./output/", outName,"/", subname,"/"), outName = "majorID_pertyName", 
       pValCutoff = 0.01, saveGeneList = T, addLabs = "", labsHide = "^ENSCAFG"
                  )


### Fig extra: upset plot of dge results by major cell types
files <- list.files(path = "/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/allCells/n3n4/linDEG", pattern=".csv", all.files=FALSE,
                        full.names=T)

df.list <- lapply(files, read.csv, header = T)

# feats.list <- lapply(df.list, function(x){feats <- x %>% filter(avg_log2FC > 0) %>% select(X)})
feats.list <- lapply(df.list, function(x){feats <- x %>% filter(avg_log2FC > 0, p_val_adj < 0.01) %>% select(X)})
feats.list <- lapply(df.list, function(x){feats <- x %>% filter(avg_log2FC < 0, p_val_adj < 0.01) %>% select(X)})


library(UpSetR)
upSet.df <- as.data.frame(unique(c(feats.list[1][[1]]$X,feats.list[2][[1]]$X,feats.list[3][[1]]$X,feats.list[4][[1]]$X,feats.list[5][[1]]$X,feats.list[5][[1]]$X,feats.list[6][[1]]$X,feats.list[7][[1]]$X)))
colnames(upSet.df) <- "gene"

upSet.df$`B cell` <- as.integer(ifelse(upSet.df$gene %in% feats.list[1][[1]]$X, 1, 0))
upSet.df$`Cycling T cell` <- as.integer(ifelse(upSet.df$gene %in% feats.list[2][[1]]$X, 1, 0))
upSet.df$Epithelial <- as.integer(ifelse(upSet.df$gene %in% feats.list[3][[1]]$X, 1, 0))
upSet.df$`Mast cell` <- as.integer(ifelse(upSet.df$gene %in% feats.list[4][[1]]$X, 1, 0))
upSet.df$Myeloid <- as.integer(ifelse(upSet.df$gene %in% feats.list[5][[1]]$X, 1, 0))
upSet.df$`Plasma cell` <- as.integer(ifelse(upSet.df$gene %in% feats.list[6][[1]]$X, 1, 0))
upSet.df$`T cell` <- as.integer(ifelse(upSet.df$gene %in% feats.list[7][[1]]$X, 1, 0))

# Plot sample distance heatmap with ComplexHeatmap
png(file = paste0("./output/", outName,"/", subname,"/",subname, "_upSet.png"), width=2200, height=2000, res=400)
par(mfcol=c(1,1))     
p <- upset(upSet.df, sets = colnames(upSet.df)[-1], cutoff = NULL,  nintersects = 7,empty.intersections = T)
p
dev.off()
