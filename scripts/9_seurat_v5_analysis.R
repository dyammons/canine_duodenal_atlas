#!/usr/bin/Rscript

#load custom functions & packages
# singularity run -B $PWD/../../../ --fakeroot ../software/seuratv5
# export TMPDIR='/scratch/alpine/dyammons@colostate.edu/cache/tmp'; python3 -m pip --cache-dir /scratch/alpine/dyammons@colostate.edu/cache/pip2 install scvi-tools
# singularity shell -B $PWD/../../../ --fakeroot ../software/shortcake.2.0.0.sif
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions_Seuratv5.R")
library(Seurat)
library(reticulate)
library(tidyverse)
library(SeuratWrappers)
library(patchwork)
library(miloR)
# library(BiocParallel)
library(SingleCellExperiment)
library(msigdbr)
library(clusterProfiler)
library(scales)

### Analysis note: 
# This script completes supplemental analysis using scVI for each subset. The T cell data is used in the manuscript
# but the rest is extra analysis not included in the paper. Generates figures associated with Supplemental figure 7.

######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Function definition   ######## <<<<<<<<<<<<<<<<<<<
######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

runSCVIintegration <- function(
    seu.obj = NULL,
    outName = ""
){
    seu.obj[["RNA"]] <- as(object = seu.obj[["RNA"]], Class = "Assay5")

    #integrate the data using all of the four Seurat v5 integration methods
    seu.obj <- integrateData(
        seu.obj = seu.obj, dout = "../output/s2/", outName = outName, 
        new.reduction.name = "integrated_v5", runAllMethods = F, 
        method = "scVIIntegration", normalization.method = "LogNormalize", 
        indReClus = T, conda_env = "/opt/conda/envs/scvi-scgen-scmomat-unitvelo"
    )

    #complete data visualization
    seu.obj <- dataVisUMAP(
        seu.obj = seu.obj, outDir = "../output/s3/", 
        outName = paste0(outName, "_integrated_v5"), final.dims = 30, 
        final.res = 0.4, stashID = "clusterID", algorithm = 3, min.dist = 0.5, 
        n.neighbors = 40, prefix = "RNA_snn_res.", assay = "RNA", 
        reduction = "integrated_v5", saveRDS = T, return_obj = T, returnFeats = T,
        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                     "CD4", "MS4A1", "PPBP","HBM")
    )
    return(seu.obj)
}

##################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Run SCVI integration - All Cells   ######## <<<<<<<<<<<<<
##################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<

#load in all cells dataset and integrate using SCVI
seu.obj.all <- readRDS("./output/s3/canine_duodenum_annotated.rds")
seu.obj <- runSCVIintegration(seu.obj = seu.obj.all, outName = "allCells_v5")

#load in SCVI integrated all cells dataset
seu.obj <- readRDS("../output/s3/allCells_v5_integrated_v5_res0.4_dims30_dist0.5_neigh40_S3.rds")
outName <- "allCells_v5"

seu.obj$celltype.l3 <- factor(ifelse(seu.obj$clusterID_2_1 == "29",
       "MS4A1_T_cell",
       as.character(seu.obj$celltype.l3)
      ))


### Fig 1a - create UMAP by major cell types
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated_v5", 
              group.by = "clusterID_integrated_v5",
              pt.size = 0.25,
              label = T,
              label.box = T
 ) + NoLegend()
p <- formatUMAP(plot = pi, smallAxes = T)
ggsave(paste("../output/", outName,  "/", outName, "_clusID.png", sep = ""), width = 7, height = 7)


### Fig 1a - create UMAP by major cell types
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated_v5", 
              group.by = "majorID",
              cols = c("#FF89B3", "#C89504", "#00ABFD", "#B983FF", "#FA7476", "#0A9B9F", "#9B8EFF"),
              pt.size = 0.25,
              label = F,
              label.box = F
 ) + NoLegend()
p <- formatUMAP(plot = pi, smallAxes = T)
ggsave(paste("../output/", outName,  "/", outName, "_fig1a.png", sep = ""), width = 7, height = 7)

### Fig supp 1a - create UMAP by celltype.l3
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated_v5", 
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
ggsave(paste("../output/", outName, "/", outName, "_sup1a.png", sep = ""), width = 7, height = 7)


### Use miloR to further validate
# Set up metadata
da_design <- as.data.frame(list(
    "Sample" = factor(c("CIE_2", "CIE_3", "CIE_4", "CIE_6", "H_1", "H_2", "H_3")),
    "Condition" = factor(c("CIE", "CIE", "CIE", "CIE", "H", "H", "H"), levels = c("H", "CIE"))
))
p <- runMilo(seu.obj = seu.obj, da_design = da_design, subName = "CIE_vs_H", blocked = F, alpha = 0.1)
p0 <- p[[1]] + ggtitle("CIE versus Healthy") + theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave(paste0("../output/", outName, "/", outName, "_milo_test.png"), width = 7, height = 7)

milo.obj <- p[[2]] 
rownames(da_design) <- da_design$Sample
da_design <- da_design[colnames(nhoodCounts(milo.obj)), , drop = FALSE]
da_results <- testNhoods(milo.obj, design = ~ Condition, design.df = da_design)
da_results %>% arrange(SpatialFDR) %>% filter(SpatialFDR < 0.4) %>% nrow()


p <- plotNhoodGraphDA(milo.obj, da_results[!is.na(da_results$logFC), ],
                          subset.nhoods = !is.na(da_results$logFC), alpha = 0.4)
p0 <- p + ggtitle("CIE versus Healthy") + theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave(paste("../output/", outName, "/", outName, "_milo_test.png", sep = ""), width = 7, height = 7)


################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Run SCVI integration - Myeloid   ######## <<<<<<<<<<<<<<
################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in processed data
seu.obj <- readRDS("./output/s3/230829_myeloid_duod_h3c4_NoIntrons_2500_res0.5_dims40_dist0.25_neigh25_S3.rds")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))
colz.df <- read.csv("./cellColz.csv", header = F)
colz.df <- colz.df[colz.df$V2 == "myeloid", ]
outName <- "myeloid"

#stash cell type IDs
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "Monocyte (c2)", "1" = "Neutrophil (c0)", 
                                   "2" = "Eosinophil (c1)", "3" = "Macrophage (c3)",
                                   "4" = "Neutrophil (c0)", "5" = "Eosinophil (c1)",
                                   "6" = "cDC1 (c4)", "7" = "IL22RA2_DC (c5)",
                                   "8" = "Eosinophil (c1)")
                       )

seu.obj$majorID_sub <- Idents(seu.obj)
seu.obj$majorID_sub <- factor(seu.obj$majorID_sub, levels = c("Neutrophil (c0)", "Eosinophil (c1)", "Monocyte (c2)", 
                                                              "Macrophage (c3)", "cDC1 (c4)", "IL22RA2_DC (c5)"))

#stash the numerical ID
clusterID_final <- table(seu.obj$majorID_sub) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
mutate(clusterID_final=row_number()-1) %>% arrange(clusterID_final) 

newID <- clusterID_final$clusterID_final
names(newID) <- clusterID_final$Var1
Idents(seu.obj) <- "majorID_sub"
seu.obj <- RenameIdents(seu.obj, newID)
table(Idents(seu.obj))
seu.obj$clusterID_final <- Idents(seu.obj)

seu.obj <- runSCVIintegration(seu.obj = seu.obj, outName = "myeloid_v5")


### Fig supp 1a - create UMAP by previous annotations
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated_v5", 
              group.by = "majorID_sub",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = TRUE,
              repel = T
)
p <- formatUMAP(plot = pi) + NoLegend() + theme(plot.title = element_text(size = 18, vjust = 1),
                                                axis.title = element_blank(),
                                                panel.border = element_blank()) + ggtitle("Myeloid canine duodenum atlas")
ggsave(paste("../output/", outName, "/", outName, "_sup1a.png", sep = ""), width = 7, height = 7)


### Use miloR to further validate
# Set up metadata
da_design <- as.data.frame(list(
    "Sample" = factor(c("CIE_1", "CIE_2", "CIE_3", "CIE_4", "H_1", "H_2", "H_3")),
    "Condition" = factor(c("CIE", "CIE", "CIE", "CIE", "H", "H", "H"), levels = c("H", "CIE"))
))
p <- runMilo(seu.obj = seu.obj, da_design = da_design, subName = "CIE_vs_H", blocked = F, alpha = 0.1)

p0 <- p[[1]] + ggtitle("CIE versus Healthy") + theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave(paste0("../output/", outName, "/", outName, "_milo_test.png"), width = 7, height = 7)

milo.obj <- p[[2]] 
rownames(da_design) <- da_design$Sample
da_design <- da_design[colnames(nhoodCounts(milo.obj)), , drop = FALSE]
da_results <- testNhoods(milo.obj, design = ~ Condition, design.df = da_design)
da_results %>% arrange(SpatialFDR) %>% filter(SpatialFDR < 0.4) %>% nrow()


p <- plotNhoodGraphDA(milo.obj, da_results[!is.na(da_results$logFC), ],
                          subset.nhoods = !is.na(da_results$logFC), alpha = 0.4)
p0 <- p + ggtitle("CIE versus Healthy") + theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave(paste("../output/", outName, "/", outName, "_milo_test.png", sep = ""), width = 7, height = 7)


################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Run SCVI integration - T cells   ######## <<<<<<<<<<<<<<
################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in processed data
seu.obj <- readRDS("./output/s3/230913_tcell_duod_h3c4_NoIntrons_2500_res0.6_dims35_dist0.3_neigh30_S3.rds")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))
colz.df <- read.csv("./cellColz.csv", header = F)
colz <- c("#C47AEA", "#EF98CB", "#FF1E6D", "#FF007F")
outName <- "tcell"

#remane idents to match the results of clustering at a resolution of 0.2 (as determined using clustree)
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

#stash the numerical ID
clusterID_final <- table(seu.obj$majorID_sub) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
mutate(clusterID_final=row_number()-1) %>% arrange(clusterID_final) 

newID <- clusterID_final$clusterID_final
names(newID) <- clusterID_final$Var1
Idents(seu.obj) <- "majorID_sub"
seu.obj <- RenameIdents(seu.obj, newID)
table(Idents(seu.obj))
seu.obj$clusterID_final <- Idents(seu.obj)

seu.obj <- runSCVIintegration(seu.obj = seu.obj, outName = "tcell_v5")

seu.obj <- readRDS("../output/s3/tcell_integrated_v5_res0.6_dims30_dist0.3_neigh30_S3.rds")
seu.obj.all <- readRDS("./output/s3/canine_duodenum_annotated.rds")
seu.obj <- AddMetaData(seu.obj, seu.obj.all$celltype.l3, col.name = "celltype.l3")
outName <- "tcell"

#correct naming change
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

Idents(seu.obj) <- "majorID_sub"
seu.obj <- RenameIdents(seu.obj, c("Non-resident" = "IL7R_high",
                                   "Tissue resident" = "GZMA_high"
                                  )
                       )
seu.obj$majorID_sub <- Idents(seu.obj)

# seu.obj <- dataVisUMAP(
#     seu.obj = seu.obj, outDir = "../output/s3/", 
#     outName = paste0(outName, "_integrated_v5"), final.dims = 30, 
#     final.res = 0.6, stashID = "clusterID_v5", algorithm = 3, min.dist = 0.3, 
#     n.neighbors = 30, prefix = "RNA_snn_res.", assay = "RNA", 
#     reduction = "integrated_v5", saveRDS = F, return_obj = T, returnFeats = T,
#     features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
#                  "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
#                  "CD4", "MS4A1", "PPBP","HBM")
# )

vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_integrated_v5",
          outName = "tell_v5", outDir = paste0("../output/violin/", outName, "/"))


### Fig supp 1a - create UMAP by previous annotations
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated_v5", 
              group.by = "clusterID_v5_integrated_v5",
              pt.size = 0.25,
              label = T,
              label.box = T
             ) + NoLegend()
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.8)  +
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5),
        axis.title = element_blank(),
        panel.border = element_blank()
    ) +
    ggtitle("scVI unsupervised clustering")
ggsave(paste0("../output/", outName, "/", outName, "_fig3d.png"), width = 7, height = 7)


pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated_v5", 
              group.by = "celltype.l3",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = TRUE,
              repel = T
) + NoLegend() + 
    labs(
        title = "T cell canine duodenum atlas",
        subtitle = "(High resolution annotations)"
    )

p <- formatUMAP(plot = pi) +
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5),
        axis.title = element_blank(),
        panel.border = element_blank()
    )

ggsave(paste("../output/", outName, "/", outName, "_sup1a.png", sep = ""), width = 7, height = 7)


pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated_v5", 
              group.by = "celltype.l3",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = TRUE,
              repel = T
) + NoLegend() + 
    labs(
        title = "T cell canine duodenum atlas"
    )

p <- formatUMAP(plot = pi, smallAxes = T) &
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        panel.border = element_blank()
    )

ggsave(paste("../output/", outName, "/", outName, "_sup1a.png", sep = ""), width = 7, height = 7)



### Fig supp 1a - create UMAP by previous annotations
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated_v5", 
              group.by = "celltype.l3",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = TRUE,
              repel = T
) + NoLegend() + 
    labs(
        title = "T cell canine duodenum atlas",
        subtitle = "(High resolution annotations)"
    )

p <- formatUMAP(plot = pi) +
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5),
        axis.title = element_blank(),
        panel.border = element_blank()
    )

ggsave(paste("../output/", outName, "/", outName, "_sup1a.png", sep = ""), width = 7, height = 7)


### Fig supp 1a - create UMAP by previous annotations
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated_v5", 
              group.by = "majorID_sub",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = TRUE,
              repel = T
) + NoLegend() + 
    labs(
        title = "T cell canine duodenum atlas",
        subtitle = "(Low resolution annotations)"
    )

p <- formatUMAP(plot = pi, smallAxes = T) &
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5)
    )
ggsave(paste("../output/", outName, "/", outName, "_sup1a.png", sep = ""), width = 7, height = 7)


### Fig supp 5c - transfer TRDC data from ROS data
seu.obj.ros <- readRDS(file = "./output/s3/230828_duod_h3_ros_res0.4_dims50_dist0.4_neigh30_S3.rds")

seu.obj.Hsub <- subset(seu.obj,
                  subset = 
                  cellSource ==  "Healthy"
                 ) 

TRDC_ex <- FetchData(seu.obj.ros, vars = "TRDC")

rownames(TRDC_ex) <- ifelse(grepl("1_1",rownames(TRDC_ex)), paste0(substr(rownames(TRDC_ex), 1, nchar(rownames(TRDC_ex))-3),"1_7") ,rownames(TRDC_ex))
rownames(TRDC_ex) <- ifelse(grepl("1_2",rownames(TRDC_ex)), paste0(substr(rownames(TRDC_ex), 1, nchar(rownames(TRDC_ex))-3),"1_8") ,rownames(TRDC_ex))
rownames(TRDC_ex) <- ifelse(grepl("1_3",rownames(TRDC_ex)), paste0(substr(rownames(TRDC_ex), 1, nchar(rownames(TRDC_ex))-3),"1_9") ,rownames(TRDC_ex))

seu.obj.Hsub <- AddMetaData(seu.obj.Hsub, metadata = TRDC_ex, col.name = "TRDC")

features = c("TRDC", "ENSCAFG00000030206")

p <- prettyFeats(seu.obj = seu.obj.Hsub, nrow = 2, ncol = 1, reduction = "umap.integrated_v5",  features = features, color = "black", order = F) 
ggsave(paste("../output/", outName, "/", outName, "_supp5c.png", sep = ""), width = 5, height = 8)


### Fig 3b - dot plot of key t cell features
features <- c("CD4", "IL7R", "GZMA","CCL4","IL17RB",
              "CD8A", "TCF7","GZMB", "GZMK", "IL13")
p <- prettyFeats(seu.obj = seu.obj, nrow = 2, ncol = 5, features = features, reduction = "umap.integrated_v5",
                 color = "black", order = F, pt.size = 0.0000001, title.size = 16, noLegend = T)
ggsave(paste("../output/", outName, "/",outName, "_fig3b.png", sep = ""), width = 12.5, height = 5)


### Use miloR to further validate
# Set up metadata
seu.obj$Sample <- seu.obj$name2
seu.obj$Condition <- factor(seu.obj$cellSource, levels = c("CIE", "H"))

da_design <- as.data.frame(list(
    "Sample" = factor(c("CIE_1", "CIE_2", "CIE_3", "CIE_4", "H_1", "H_2", "H_3")),
    "Condition" = factor(c("CIE", "CIE", "CIE", "CIE", "H", "H", "H"), levels = c("H", "CIE"))
))

miloOut <- runMilo(
    seu.obj = seu.obj, da_design = da_design, groupBy = "Sample",
    splitBy = "Condition", outName = outName, subName = "CIE_vs_H",
    embedding = seu.obj@reductions$umap.integrated_v5@cell.embeddings,
    blocked = FALSE, alpha = 0.2
)


p1 <- miloOut[[1]] + 
    ggtitle("CIE versus healthy") + 
    theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    ) + 
    scale_fill_gradient2(
        name = "log2(FC)",
        low = "blue", mid = "white", high = "red", midpoint = 0
    ) + 
    scale_edge_width(
        name = "Overlap size", range = c(0.2, 3)
    ) + 
    guides(
        size = guide_legend(title = "Neighborhood\nsize",
                            override.aes = list(fill = NA, shape = 21))
    )
ggsave(paste0("../output/", outName, "/", outName, "_milo_sig.png"), width = 8, height = 7)


################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Run SCVI integration - Epithelial   ######## <<<<<<<<<<<<<<
################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in preprocessed data
seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/s3/230916_duod_duod_h3c4_NoIntrons_2500_res1.7_dims40_dist0.25_neigh25_S3.rds")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))
colz.df <- read.csv("./cellColz.csv", header = F)
colz.df <- colz.df[colz.df$V2 == "duod", ]
outName <- "epithelial"

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

seu.obj <- runSCVIintegration(seu.obj = seu.obj, outName = "epithelial_v5")

### reload here
seu.obj.all <- readRDS("./output/s3/canine_duodenum_annotated.rds")
seu.obj <- readRDS("../output/s3/epithelial_v5_integrated_v5_res0.4_dims30_dist0.5_neigh40_S3.rds")
seu.obj <- AddMetaData(seu.obj, seu.obj.all$celltype.l3, col.name = "celltype.l3")
outName <- "epithelial"

vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_integrated_v5",
          outName = "tell_v5", outDir = paste0("../output/violin/", outName, "/"))

### Fig supp 1a - create UMAP by previous annotations
pi <- DimPlot(seu.obj, 
              reduction = "umap.integrated_v5", 
              group.by = "celltype.l3",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = TRUE,
              repel = T
)
p <- formatUMAP(plot = pi) + NoLegend() + theme(plot.title = element_text(size = 18, vjust = 1),
                                                axis.title = element_blank(),
                                                panel.border = element_blank()) + ggtitle("Epithelial canine duodenum atlas")
ggsave(paste("../output/", outName, "/", outName, "_sup1a.png", sep = ""), width = 7, height = 7)


### Use miloR to further validate
# Set up metadata
da_design <- as.data.frame(list(
    "Sample" = factor(c("CIE_1", "CIE_2", "CIE_3", "CIE_4", "H_1", "H_2", "H_3")),
    "Condition" = factor(c("CIE", "CIE", "CIE", "CIE", "H", "H", "H"), levels = c("H", "CIE"))
))
p <- runMilo(seu.obj = seu.obj, da_design = da_design, subName = "CIE_vs_H", blocked = F, alpha = 0.1)

p0 <- p[[1]] + ggtitle("CIE versus Healthy") + theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave(paste0("../output/", outName, "/", outName, "_milo_test.png"), width = 7, height = 7)

milo.obj <- p[[2]] 
rownames(da_design) <- da_design$Sample
da_design <- da_design[colnames(nhoodCounts(milo.obj)), , drop = FALSE]
da_results <- testNhoods(milo.obj, design = ~ Condition, design.df = da_design)
da_results %>% arrange(SpatialFDR) %>% filter(SpatialFDR < 0.1) %>% nrow()


p <- plotNhoodGraphDA(milo.obj, da_results[!is.na(da_results$logFC), ],
                          subset.nhoods = !is.na(da_results$logFC), alpha = 0.15)
p0 <- p + ggtitle("CIE versus Healthy") + theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
ggsave(paste("../output/", outName, "/", outName, "_milo_test.png", sep = ""), width = 7, height = 7)


################################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   End SCVI integration analysis   ######## <<<<<<<<<<<<<<
################################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<



