#!/usr/bin/Rscript

#load custom functions & packages
# source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")
library(circlize)

### Analysis note: 
# This script loads in the previously processed Seurat object (./output/s3/230816_duod_h3c4_NoIntrons_res1.3_dims40_dist0.3_neigh50_S3.rds)
# then loads all of the major cell type datasets annoated though independent reclustering and transfers the cell type labels to the full
# dataset. The script generates all figures assocaited with Figure 1.

################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Transfer cell type annotations   ######## <<<<<<<<<<<<<<
################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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

#stash the idents
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

#load in processed tcell subset data -- to get ILC2s
seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/s3/230913_tcell_duod_h3c4_NoIntrons_2500_res0.6_dims35_dist0.3_neigh30_S3.rds")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))
colz.df <- read.csv("./cellColz.csv", header = F)
outName <- "tcell"
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "Tissue resident", "1" = "Tissue resident", 
                                   "2" = "Non-resident", "3" = "Memory",
                                   "4" = "Tissue resident", "5" = "Non-resident",
                                   "6" = "Non-resident", "7" = "Non-resident",
                                   "8" = "Memory", "9" = "Tissue resident",
                                   "10" = "Tissue resident", "11" = "Non-resident",
                                   "12" = "Tissue resident", "13" = "ILC2",
                                   "14" = "Non-resident","15" = "Non-resident")
                       )
seu.obj$majorID_sub <- Idents(seu.obj)
seu.obj <- subset(seu.obj, subset = majorID_sub == "ILC2")
seu.obj$majorID_sub <- droplevels(Idents(seu.obj))
ct.l3 <- c(ct.l3,seu.obj$majorID_sub)

#load in processed tcell subset data -- to get other T cell pops
seu.obj <- readRDS("./output/s3/240112_tcell_noILC_duod_h3c4_NoIntrons_2500_res0.6_dims30_dist0.1_neigh10_S3.rds")

#clusterID_sub 10 is artifactually split across the UMAP space. Use UMAP coordinates to split out cluster
#extract the barcodes of cells below -2.5 on the y-axis of the UMAP
lowerC10 <- seu.obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% filter(rownames(.) %in% WhichCells(seu.obj, idents = "10"), UMAP_2 < -2.5)
seu.obj$clusterID_sub <- as.factor(ifelse(seu.obj$clusterID_sub == "10",
                                          ifelse(names(seu.obj$clusterID_sub) %in% rownames(lowerC10), "10.2", "10.1"), as.character(seu.obj$clusterID_sub)
                                         )
                                  )
table(seu.obj$clusterID_sub)
#stash new cell idents
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "CD8_eff", "1" = "CD8_eff", 
                                   "2" = "Tnaive", "3" = "CD8_TRM",
                                   "4" = "gdT_1", "5" = "Tnaive",
                                   "6" = "NK_T", "7" = "Tnaive",
                                   "8" = "CD8_mem", "9" = "gdT_2",
                                   "10.1" = "Tnaive", "10.2" = "CD8_eff",
                                   "11" = "CD8_eff", "12" = "Treg",
                                   "13" = "NK","14" = "T_IFN")
                       )
seu.obj$majorID_sub_noILC <- Idents(seu.obj)
ct.l3 <- c(ct.l3,seu.obj$majorID_sub_noILC)

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
seu.obj$celltype.l3 <- as.factor(ifelse(is.na(seu.obj$celltype.l3), "remove", as.character(seu.obj$celltype.l3)))
seu.obj <- subset(seu.obj, invert = T,
                  subset = celltype.l3 == "remove"
                 )
seu.obj$celltype.l3 <- droplevels(seu.obj$celltype.l3)
saveRDS(seu.obj, "./output/s3/canine_duodenum_annotated.rds")


########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin all cells analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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


pi <- DimPlot(seu.obj,
        reduction = "umap", 
        cols = "grey",
        pt.size = 0.5,
        label = F,
        cells.highlight = WhichCells(seu.obj, expression = nCount_RNA < 1000),
        cols.highlight = "orchid",
        order = F
        #label.box = TRUE
 ) + NoLegend()
p <- formatUMAP(pi) 
ggsave(paste("./output/lowQ_highlight.png", sep = ""), width =7, height = 7)

filter_inc <- rep("exclude", length(WhichCells(seu.obj, expression = nCount_RNA < 1000)))
names(filter_inc) <- WhichCells(seu.obj, expression = nCount_RNA < 1000)
seu.obj <- AddMetaData(seu.obj, filter_inc, col.name = "exclude")
df <- table(seu.obj$exclude, seu.obj$celltype.l3) %>%
    as.data.frame() %>%
    left_join(as.data.frame(table(seu.obj$celltype.l3)), by = c("Var2" = "Var1")) %>%
    mutate(
        pctRemoved = round(( Freq.x / Freq.y ) * 100, 2)
    ) %>%
    arrange(desc(pctRemoved)) %>%
    select(-Var1, -Freq.x, -Freq.y) %>%
    rename(`Cell type` = Var2, `Additional % removed` = pctRemoved)


### Fig extra - check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


### Supp data - run FindAllMarkers on major and high-res cell types

#run on major ID
saveName <- "240128_allCells_duod_h3c4_NoIntrons_2500"
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_pertyName", numOfFeats = 24, outName = saveName,
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), returnViln = F
                    )

#export surface marker data with finaallmarkers
surface.markers <- read.csv("./surface_master.csv")[ ,c("UniProt.gene", "UniProt.description", "Surfaceome.Label", "Surfaceome.Label.Source")] %>% filter(!duplicated(UniProt.gene))
cluster.markers <- read.csv("./output/viln/allCells/240128_allCells_duod_h3c4_NoIntrons_2500_gene_list.csv")
write.csv(cluster.markers[ ,c(7,8,2:6)] %>% left_join(surface.markers, by = c("gene" = "UniProt.gene")),
          file = "./output/supplementalData/supplemental_data_1.csv", row.names = F)

#run on high-res dataset
saveName <- "240128_allCells_duod_h3c4_NoIntrons_2500"
vilnPlots(seu.obj = seu.obj, groupBy = "celltype.l3", numOfFeats = 24, outName = saveName,
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), returnViln = F
                    )

#export surface marker data with finaallmarkers
surface.markers <- read.csv("./surface_master.csv")[ ,c("UniProt.gene", "UniProt.description", "Surfaceome.Label", "Surfaceome.Label.Source")] %>% filter(!duplicated(UniProt.gene))
cluster.markers <- read.csv("./output/viln/allCells/240128_allCells_duod_h3c4_NoIntrons_2500_gene_list.csv")
write.csv(cluster.markers[ ,c(7,8,2:6)] %>% left_join(surface.markers, by = c("gene" = "UniProt.gene")),
          file = "./output/supplementalData/supplemental_data_2.csv", row.names = F)


#export data for cellbrowser
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
              label = T,
              label.box = T,
              shuffle = TRUE,
              repel = T
)
p <- formatUMAP(plot = pi) + NoLegend() + theme(plot.title = element_text(size = 18, vjust = 1),
                                                axis.title = element_blank(),
                                                panel.border = element_blank()) + ggtitle("Annotated canine duodenum atlas")
ggsave(paste("./output/", outName, "/", outName, "_sup1a.png", sep = ""), width = 7, height = 7)

### Fig 1a: create labels to be cropped onto UMAP
majorColors.df <- as.data.frame(c("T cell", "Epithelial", "Myeloid", "Plasma cell", "Cycling T cell", "Mast", "B cell"))
colnames(majorColors.df) <- "clusterID_final"
majorColors.df$colz <- c("#FF89B3", "#C89504", "#00ABFD", "#B983FF", "#FA7476", "#0A9B9F", "#9B8EFF")
majorColors.df$title <- "All cells"
majorColors.df$labCol <- "black"
majorColors.df$labz <- 0:6
leg <- cusLeg(legend = majorColors.df, clusLabel = "labz", legLabel = "clusterID_final", colorz = "colz",labCol = "labCol",colz = 1, rowz = NULL, groupLabel = "title", dotSize = 10, groupBy = "title",sortBy = "labz", compress_x = 0.9, textSize = 6)
ggsave(paste("./output/", outName, "/", outName, "_leg_forUMAP_labels.png", sep = ""), width = 2.75, height = 7)


### Fig 1a - create UMAP by major cell types
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_final",
              cols = c("#FF89B3", "#C89504", "#00ABFD", "#B983FF", "#FA7476", "#0A9B9F", "#9B8EFF"),
              pt.size = 0.25,
              label = F,
              label.box = F
 ) + NoLegend()
p <- formatUMAP(plot = pi, smallAxes = T)
ggsave(paste("./output/", outName,  "/", subname,"/", outName, "_fig1a.png", sep = ""), width = 7, height = 7)


### Fig supp 1b- reference map using PBMC data
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


### Fig supp 2c - reference map using human epitheial data
reference <- MuDataSeurat::ReadH5AD("./epi_log_counts02_v2.h5ad")

#prep ref
reference <- SCTransform(reference, verbose = FALSE)
reference <- RunPCA(reference)

#map data
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

### Fig supp 1d - reference map using human t cell data
reference <- MuDataSeurat::ReadH5AD("./Tcell_log_counts02_v2.h5ad")

#prep ref
reference <- SCTransform(reference, verbose = FALSE)
reference <- RunPCA(reference)

anchors <- FindTransferAnchors(
    reference = reference,
    query = seu.obj,
    normalization.method = "SCT",
    reference.reduction = "pca", 
    dims= 1:50
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
                                                        )  + guides(fill = guide_legend(nrow = 2))
ggsave(paste("./output/", outName,"/", subname, "/",outName, "_fig1d.png", sep = ""), width = 5, height = 4)


### Fig supp 2a: umap by sample
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
ggsave(paste("./output/", outName, "/", subname, "/", outName, "_sup2a.png", sep = ""), width =7, height = 7)


### Fig supp 2b: stats by cie vs healthy
freqy <- freqPlots(seu.obj, method = 2, nrow= 2, groupBy = "majorID_pertyName", legTitle = "Cell source",refVal = "name2", showPval = T,
              namez = "name2", 
              colz = "colz"
              )
ggsave(paste("./output/", outName, "/", subname, "/",outName, "_supp2b.png", sep = ""), width = 8, height = 6)

#run levene test to ensure proper test was used
res.ftest <- lapply(levels(freqy$data$majorID_pertyName), function(x){
    data.df <- freqy$data[freqy$data$majorID_pertyName ==  x, ]
    car::leveneTest(freq ~ cellSource, data = data.df)
})
res.ftest
               


### Use miloR to further validate
library(miloR)
library(BiocParallel)

runMilo <- function(
    seu.obj = NULL, 
    da_design = NULL, 
    subName = "", 
    blocked = TRUE, 
    ...
    ){
    
    # Hard coded metadata re-naming
    seu.obj$Sample <- seu.obj$name2
    seu.obj$Condition <- seu.obj$cellSource

    # Convert from Seurat to sce object
    sce <- as.SingleCellExperiment(seu.obj)
    reducedDim(sce, "PCA") <- seu.obj@reductions$pca@cell.embeddings
    reducedDim(sce, "UMAP") <- seu.obj@reductions$umap@cell.embeddings

    # Preprocess using miloR to ID neighboorhoods
    milo.obj <- Milo(sce)
    milo.obj$Sample <- droplevels(factor(milo.obj$Sample))
    milo.obj <- buildGraph(milo.obj, k = 30, d = 40)
    milo.obj <- makeNhoods(milo.obj, prop = 0.2, k = 30, d = 40, refined = TRUE, refinement_scheme = "graph")
    p <- plotNhoodSizeHist(milo.obj)
    ggsave(paste0("./output/", outName, "/", outName, "_NhoodSize.png"), width = 7, height = 7)
    
    milo.obj <- countCells(milo.obj, meta.data = data.frame(colData(milo.obj)), samples = "Sample")

    # Set up metadata
    rownames(da_design) <- da_design$Sample
    da_design <- da_design[colnames(nhoodCounts(milo.obj)), , drop = FALSE]

    # Calc distance between neighborhoods and test for DA
    milo.obj <- calcNhoodDistance(milo.obj, d = 40)
    if(blocked){
        da_results <- testNhoods(milo.obj, design = ~ Dog + Condition, design.df = da_design)
    } else{
        da_results <- testNhoods(milo.obj, design = ~ Condition, design.df = da_design)
    }
    
    n_diff <- da_results %>% arrange(SpatialFDR) %>% filter(SpatialFDR < 0.1) %>% nrow()
    if(n_diff == 0){
        message(
            paste(
                "No differentially abundant Nhoods at alpha = 0.1!",
                "The lowest spaitally adjusted P value is:", min(da_results$SpatialFDR),
                "\n Although not reccomended, you can increase alpha to the lowest",
                "spaitally adjusted P value to get an idea of which regoins of the",
                "UMAP are trendy."
            )
        )
    }

    # Plot the results (by neighborhood)
    milo.obj <- buildNhoodGraph(milo.obj)
    p <- plotNhoodGraphDA(milo.obj, da_results[!is.na(da_results$logFC), ],
                          subset.nhoods=!is.na(da_results$logFC), ...)
    ggsave(paste("./output/", outName, "/", subName, "_milo.png", sep = ""), width = 6, height = 6)
    return(list(p, milo.obj))
}

# Set up metadata
da_design <- as.data.frame(list(
    "Sample" = factor(c("CIE_1", "CIE_2", "CIE_3", "CIE_4", "H_1", "H_2", "H_3")),
    "Condition" = factor(c("CIE", "CIE", "CIE", "CIE", "H", "H", "H"), levels = c("H", "CIE"))
))
p <- runMilo(seu.obj = seu.obj, da_design = da_design, subName = "CIE_vs_H", blocked = F, alpha = 0.1)
p0 <- p[[1]] + ggtitle("CIE versus Healthy") + theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave(paste0("./output/", outName, "/", outName, "_milo_test.png"), width = 7, height = 7)

milo.obj <- p[[2]] 
rownames(da_design) <- da_design$Sample
da_design <- da_design[colnames(nhoodCounts(milo.obj)), , drop = FALSE]
da_results <- testNhoods(milo.obj, design = ~ Condition, design.df = da_design)
da_results %>% arrange(SpatialFDR) %>% filter(SpatialFDR < 0.9) %>% nrow()


p0 <- plotNhoodGraphDA(milo.obj, da_results[!is.na(da_results$logFC), ],
                          subset.nhoods = !is.na(da_results$logFC), alpha = 1)
p0 <- p0 + ggtitle("CIE versus Healthy") + theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave(paste("./output/", outName, "/", outName, "_milo_test.png", sep = ""), width = 7, height = 7)




### Fig 1e: dge analysis all cells
seu.obj$allCells <- "DGE analysis of all cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
linDEG(seu.obj = seu.obj, threshold = 1, thresLine = F, groupBy = "allCells", comparision = "cellSource", contrast = c("CIE", "Healthy"),
       outDir = paste0("./output/", outName,"/", subname,"/"), 
       outName = "all_cells", cluster = NULL, labCutoff = 10, noTitle = F, pValCutoff = 0.01, logfc.threshold = 0.58,
                   colUp = "red", colDwn = "blue", subtitle = T, returnUpList = F, returnDwnList = F, forceReturn = F, useLineThreshold = F, saveGeneList = T, addLabs = ""
                  )


### Fig extra: dge analysis by major cell types
linDEG(seu.obj = seu.obj, groupBy = "majorID_pertyName", comparision = "cellSource", contrast= c("CIE","Healthy"),
       outDir = paste0("./output/", outName,"/", subname,"/linDEG/"), outName = "majorID_pertyName", 
       pValCutoff = 0.01, logfc.threshold = 0.58, saveGeneList = T, addLabs = "",
                  )


### Fig extra - pseudobulk DEGs btwn healthy and cie
seu.obj$allCells <- "All cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
createPB(seu.obj = seu.obj, groupBy = "allCells", comp = "cellSource", biologicalRep = "name2", lowFilter = T, 
         dwnSam = F, min.cell = 15, clusters = NULL, outDir = paste0("./output/", outName,"/pseudoBulk/") , 
         grepTerm = "H", grepLabel = c("Healthy","CIE") #improve - fix this so it is more functional
        )

p <- pseudoDEG(metaPWD = paste0("./output/", outName,"/pseudoBulk/allCells_deg_metaData.csv"), returnDDS = F, 
               padj_cutoff = 0.1, lfcCut = 1, outDir = paste0("./output/", outName,"/pseudoBulk/"), 
               outName = "allCells", idents.1_NAME = "CIE", idents.2_NAME = "Healthy", 
               inDir = paste0("./output/", outName,"/pseudoBulk/"), title = "All cells", fromFile = T, meta = NULL, 
               pbj = NULL, returnVolc = T, paired = F, pairBy = "", minimalOuts = F, saveSigRes = T, 
               filterTerm = "^ENSCAF", addLabs = NULL, mkDir = T, strict_lfc = F, labSize = 4.5
              )
pi <- p[[1]] + scale_y_continuous(limits = c(0, 3))
pi  <- prettyVolc(plot = pi, rightLab = "Up in CIE", leftLab = "Up in healthy", arrowz = T, lfcCut = 1, y = 2.7) + labs(title = "CIE vs Healthy (within all cells)", x = "log2(Fold change)") + NoLegend() + theme(panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 2),
                                      axis.line = element_blank(),
                                      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 2))
ggsave(paste("./output/", outName, "/", outName, "_supp6e.png", sep = ""), width = 7, height = 7)


### Fig extra - pseudobulk DEGs btwn healthy and cie
createPB(seu.obj = seu.obj, groupBy = "celltype.l3", comp = "cellSource", biologicalRep = "name2", lowFilter = T, 
         dwnSam = F, min.cell = 15, clusters = NULL, outDir = paste0("./output/", outName,"/pseudoBulk/") , 
         grepTerm = "H", grepLabel = c("Healthy","CIE") #improve - fix this so it is more functional
        )

p <- pseudoDEG(metaPWD = paste0("./output/", outName,"/pseudoBulk/celltype.l3_deg_metaData.csv"), returnDDS = F, 
               padj_cutoff = 0.1, lfcCut = 1, outDir = paste0("./output/", outName,"/pseudoBulk/"), 
               outName = "allCells", idents.1_NAME = "CIE", idents.2_NAME = "Healthy", 
               inDir = paste0("./output/", outName,"/pseudoBulk/"), title = "All cells", fromFile = T, meta = NULL, 
               pbj = NULL, returnVolc = T, paired = F, pairBy = "", minimalOuts = F, saveSigRes = T, 
               filterTerm = "^ENSCAF", addLabs = NULL, mkDir = T, strict_lfc = F, labSize = 4.5
              )


### Supp fig xx -- heatmap of degs by each cluster
files <- lapply(levels(seu.obj$celltype.l3), function(x){
    list.files(path = paste0("./output/allCells/pseudoBulk/", x), pattern = ".csv", all.files = FALSE, full.names = T)
})

df.list <- lapply(unlist(files), read.csv, header = T)

seu.obj$type <- factor(paste0(as.character(seu.obj$celltype.l3), "--", as.character(seu.obj$cellSource)),
                       levels = paste0(rep(levels(seu.obj$celltype.l3), each = 2), "--", c("Healthy", "CIE")))

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
clus <- gg_color_hue(length(levels(seu.obj$celltype.l3)))
names(clus) <- levels(seu.obj$celltype.l3)
cond_colz <- c("mediumseagreen","mediumpurple1")
names(cond_colz) <- c("Healthy","CIE")

# heat_col <- viridis(option = "magma",100)
ha <- HeatmapAnnotation(
    Cluster = factor(unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"--")[[1]][1]})),
                     levels = levels(seu.obj$celltype.l3)),
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
    row_names_gp = gpar(fontsize = 8),
    show_row_names = T,
    cluster_columns = F,
    top_annotation = ha,
    show_column_names = F,
    column_split = factor(unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"--")[[1]][1]})),
                          levels = levels(seu.obj$celltype.l3)),
    row_title = NULL,
    column_title = NULL,
    heatmap_legend_param = list(
        title = "Scaled expression",
        direction = "horizontal"
        ),
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sig.mat[i, j], x, y, gp = gpar(fontsize = 8, col = "black"))
    }
)

png(file = paste0("./output/", outName, "/", outName, "_fig3e.png"), width=3750, height=4500, res=400)
par(mfcol=c(1,1))   
draw(ht, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "bottom")

# for(i in 1:length(levels(seu.obj$majorID_sub_inc))){
#     decorate_annotation("Cluster", slice = i, {
#         grid.text(paste0("c", (1:12) - 1)[i], just = "center")
#     })
# }
dev.off()

### Fig 1b - dot plot by major cell types
seu.obj$majorID_pertyName_split <- paste0(as.character(seu.obj$majorID_pertyName), " (", as.character(seu.obj$cellSource), ")")
fig1c <- majorDot(seu.obj = seu.obj, groupBy = "majorID_pertyName_split",
                  yAxis = NULL, #c("Epithelial","T cell","Myeloid","Plasma cell","B cell","Mast cell","Cycling T cell"),
                  features = c("KCNK16", "CYP1B1", "IL17F")
                 ) + theme(axis.title = element_blank(),
                           axis.text = element_text(size = 12)) +  scale_colour_viridis(option="magma", name='Average\nexpression', 
                                                                                        breaks = c(-0.5, 1, 2),
                                                                                        labels = c("-0.5", "1", "2")
                                                                                       ) + guides(color = guide_colorbar(title = 'Scaled\nExpression  ')) + coord_flip()
ggsave(paste("./output/", outName, "/", subname, "/", outName, "_fig1b.png", sep = ""), width = 6, height = 3)


### Fig 2b - create violin plots for key feats
features <- c("KCNK16", "CYP1B1", "IL17F")

pi <- VlnPlot(object = seu.obj,
              pt.size = 0,
              same.y.lims = F,
              group.by = "majorID_pertyName",
              split.by = "cellSource",
              combine = T,
              cols = c("#FF89B3", "#C89504", "#00ABFD", "#B983FF", "#FA7476", "#0A9B9F", "#9B8EFF"),
              stack = T,
              fill.by = "ident",
              flip = T,
              features = features
             ) + NoLegend() + theme(axis.ticks = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.title.x = element_blank(),
                                    plot.margin = unit(c(7, 7, 7, 21), "pt"))
ggsave(paste("./output/", outName, "/", outName, "_fig2b.png", sep = ""), width = 5, height =6)



### Fig extra - pseudobulk DEGs btwn healthy and cie
createPB(seu.obj = seu.obj, groupBy = "majorID", comp = "cellSource", biologicalRep = "name2", lowFilter = T, 
         dwnSam = F, min.cell = 15, clusters = NULL, outDir = paste0("./output/", outName,"/pseudoBulk/") , 
         grepTerm = "H", grepLabel = c("Healthy","CIE") #improve - fix this so it is more functional
        )
p <- pseudoDEG(metaPWD = paste0("./output/", outName,"/pseudoBulk/majorID_deg_metaData.csv"), returnDDS = F, 
               padj_cutoff = 0.1, lfcCut = 1, outDir = paste0("./output/", outName,"/pseudoBulk/"), 
               outName = "allCells", idents.1_NAME = "CIE", idents.2_NAME = "Healthy", 
               inDir = paste0("./output/", outName,"/pseudoBulk/"), title = "All cells", fromFile = T, meta = NULL, 
               pbj = NULL, returnVolc = T, paired = F, pairBy = "", minimalOuts = F, saveSigRes = T, 
               filterTerm = "^ENSCAF", addLabs = NULL, mkDir = T, strict_lfc = F
              )


### Fig 1f: heatmap of dge results by major cell types
files <- lapply(levels(seu.obj$majorID), function(x){
    list.files(path = paste0("./output/allCells/pseudoBulk/", x), pattern = ".csv", all.files = FALSE, full.names = T)
})

files <- unlist(files)[-7]
df.list <- lapply(files, read.csv, header = T)

cnts_mat <- do.call(rbind, df.list)  %>% 
    mutate(
        direction = ifelse(log2FoldChange > 0, "Up", "Down")
    ) %>% 
    group_by(gs_base, direction) %>% 
    summarize(nRow = n()) %>% 
    pivot_wider(names_from = gs_base, values_from = nRow) %>% 
    as.matrix() %>% t()

colnames(cnts_mat) <- cnts_mat[1,]
cnts_mat <- cnts_mat[-c(1),]
class(cnts_mat) <- "numeric"

#order by number of total # of DEGs
orderList <- rev(rownames(cnts_mat)[order(rowSums(cnts_mat))])
cnts_mat <- cnts_mat[match(orderList, rownames(cnts_mat)),]        
cnts_mat[is.na(cnts_mat)] <- 0
rownames(cnts_mat) <- c("Myeloid", "Mast cell", "T cell", "Epithelial", "Plasma cell", "Cycling T cell")

png(file = paste0("./output/", outName, "/", subname, "/",outName, "_fig1f.png"), width=1500, height=2000, res=400)
par(mfcol=c(1,1))         
ht <- Heatmap(cnts_mat,#name = "mat", #col = col_fun,
              name = "# of DEGs",
              cluster_rows = F,
              row_title = "Cell type",
              col = circlize::colorRamp2(c(0,max(cnts_mat)), colors = c("white","red")),
              cluster_columns = F,
              column_title = gt_render(
                  paste0("<span style='font-size:18pt; color:black'># of DEGs</span><br>",
                         "<span style='font-size:12pt; color:black'>(CIE vs Healthy)</span>")
              ),
              show_column_names = TRUE,
              column_title_side = "top",
              column_names_rot = 0,
              column_names_centered = TRUE,
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                                          labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm")),
              cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(sprintf("%.0f", as.matrix(cnts_mat)[i, j]), x, y, gp = gpar(fontsize = 14, col = "black"))
              })
draw(ht, padding = unit(c(2, 12, 2, 5), "mm"),show_heatmap_legend = FALSE)
dev.off()


### Fig sup 2c: gsea of the DGE results using majorID_sub_big
df <- read.csv("./output/allCells/n3n4/linDEG/majorID_pertyName_Mast_cell_geneList.csv")
upGenes <- df %>% filter(avg_log2FC > 0) %>% pull(X)
dwnGenes <- df %>% filter(avg_log2FC < 0) %>% pull(X)
p <- plotGSEA(geneList = upGenes, geneListDwn = dwnGenes, category = "C5", subcategory = "GO:BP", 
              upCol = "blue", dwnCol = "red", size = 3)

minVal <- -5
maxVal <- 7
pi <- p + scale_x_continuous(limits = c(minVal, maxVal), name = "Signed log10(padj)") + 
    theme(axis.title=element_text(size = 16),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
         ) + 
    annotate("segment", x = -0.1, 
             y = 17, 
             xend = minVal, 
             yend = 17, 
             lineend = "round", linejoin = "bevel", linetype ="solid", colour = "blue",
             size = 1, arrow = arrow(length = unit(0.1, "inches"))
            ) + 
    annotate(geom = "text", x = (minVal-0.1*1.5)/2-0.1*1.5, 
             y = 18,
             label = "Repressed in CIE",
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
             label = "Induced in CIE",
             hjust = 0.5,
             vjust = 1.5,
             size = 5) + 
    ggtitle("GSEA within Mast cells")
ggsave(paste("./output/", outName, "/", subname, "/", outName, "_sup2c.png", sep = ""), width = 10, height = 7)


########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end all cells analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<