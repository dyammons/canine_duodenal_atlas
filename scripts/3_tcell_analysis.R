#!/usr/bin/Rscript

#load custom functions & packages
# source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")
library(circlize)

### Analysis note: 
# This script loads in the previously processed Seurat object (./output/s3/230816_duod_h3c4_NoIntrons_res1.3_dims40_dist0.3_neigh50_S3.rds)
# then subsets on T cells and generates all figures assocaited with Figure 4, Supplemental figure 6 and 8

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin T cell preprocessing   ######## <<<<<<<<<<<<<<
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


#another susspected low quality cluster identified -- remove
seu.obj.sub <- subset(seu.obj, invert = T,
                  subset = 
                  clusterID_sub ==  "13")

table(seu.obj.sub$clusterID_sub)
table(seu.obj.sub$clusterID_2_1)
table(seu.obj.sub$orig.ident)

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "./output/s2/", subName = "240112_tcell_noILC_duod_h3c4_NoIntrons_2500", 
                     preSub = T, nfeatures = 2500, vars.to.regress = "percent.mt", saveRDS = F
                    )

#clustree to determine clus resolution
# seu.obj <- readRDS(file = "./output/s2/230913_tcell_duod_h3c4_NoIntrons_2500_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "240112_tcell_noILC_duod_h3c4_NoIntrons_2500", test_dims = c("40","35", "30"), algorithm = 3, prefix = "integrated_snn_res.")

#visulize the data & evaluate
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "240112_tcell_noILC_duod_h3c4_NoIntrons_2500", final.dims = 35, final.res = 0.6, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin T cell analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in processed data
seu.obj <- readRDS("./output/s3/230913_tcell_duod_h3c4_NoIntrons_2500_res0.6_dims35_dist0.3_neigh30_S3.rds")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))
colz.df <- read.csv("./cellColz.csv", header = F)
colz <- c("#C47AEA", "#EF98CB", "#603FEF", "#FF007F")
outName <- "tcell"

#remane idents to match the results of clustering at a resolution of 0.2 (as determined using clustree)
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "GZMA_high", "1" = "GZMA_high", 
                                   "2" = "IL7R_high", "3" = "Memory",
                                   "4" = "GZMA_high", "5" = "IL7R_high",
                                   "6" = "IL7R_high", "7" = "IL7R_high",
                                   "8" = "Memory", "9" = "GZMA_high",
                                   "10" = "GZMA_high", "11" = "IL7R_high",
                                   "12" = "GZMA_high", "13" = "ILC2",
                                   "14" = "IL7R_high","15" = "IL7R_high")
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


### Fig extra - check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/",outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


### Fig supp 6d - run singleR
pp <- singleR(seu.obj = seu.obj, outName = "230913_tcell_duod_h3c4_NoIntrons_2500", clusters = "clusterID_sub", outDir = "./output/singleR/")
namez <- pp[[1]]$SingleR.human_ref1
names(namez) <- rownames(pp[[1]])

seu.obj <- AddMetaData(seu.obj, namez, col.name = "SingleR_ref1")

pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "SingleR_ref1",
              #cols = levels(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F
)
pi <- formatUMAP(plot = pi) + NoLegend() + theme(plot.title = element_text(size = 18, vjust = 1),
                                                 axis.title = element_blank(),
                                                 panel.border = element_blank()) + ggtitle("SingleR Human Primary Cell Atlas reference")
ggsave(paste0("./output/", outName, "/", outName, "_supp6d.png"), width = 7, height = 7)


### Supp fig 6b - reference map using PBMC data
reference <- readRDS(file = "../../k9_PBMC_scRNA/analysis/output/s3/final_dataSet_HvO.rds") 
#download reference with utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_final_dataSet_HvO.rds.gz", dest = "/pwd/to/dir/final_dataSet_HvO.rds.gz")

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
              shuffle = F,
              repel = T
) + NoLegend() + theme(plot.title = element_text(size = 18, hjust = 0.5)) + ggtitle("Canine PBMC reference mapping")

p <- formatUMAP(plot = pi, smallAxes = T) & theme(plot.title = element_text(size = 18, hjust = 0.5))
ggsave(paste("./output/", outName,"/",outName, "_supp6b.png", sep = ""), width = 7, height = 7)


### Supp fig 6c - reference map using gut atlas data
reference <- MuDataSeurat::ReadH5AD("./Tcell_log_counts02_v2.h5ad")
#download reference from https://www.gutcellatlas.org/#datasets

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
seu.obj$predicted.human <- seu.obj$predicted.id
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "predicted.id",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F,
              repel = T
)
pi <- formatUMAP(plot = pi) + NoLegend() + theme(plot.title = element_text(size = 18, vjust = 1),
                                                 axis.title = element_blank(),
                                                 panel.border = element_blank()) + ggtitle("Human T cell reference mapping")
ggsave(paste("./output/", outName,"/",outName, "_supp6c.png", sep = ""), width = 7, height = 7)


### supp data 5 - cell type gene signatures
#generate violin plots for each cluster
saveName <- "240128_tcell_duod_h3c4_NoIntrons_2500"
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_sub", numOfFeats = 24, outName = saveName,
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )

#export surface marker data with FindAllMarkers
surface.markers <- read.csv("./surface_master.csv")[ ,c("UniProt.gene", "UniProt.description", "Surfaceome.Label", "Surfaceome.Label.Source")] %>% filter(!duplicated(UniProt.gene))
cluster.markers <- read.csv("./output/viln/tcell/240128_tcell_duod_h3c4_NoIntrons_2500_gene_list.csv")
write.csv(cluster.markers[ ,c(7,8,2:6)] %>% left_join(surface.markers, by = c("gene" = "UniProt.gene")),
          file = "./output/supplementalData/supplemental_data_5.csv", row.names = F)


### Fig 4a - tcell colorized by annotation clustering
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_final",
              pt.size = 0.25,
              cols = colz,
              label = F,
              label.box = F
             ) + NoLegend()
p <- formatUMAP(plot = pi, smallAxes = T)
ggsave(paste0("./output/", outName, "/", outName, "_fig4a.png"), width = 7, height = 7)


### Fig 4b - feature plot of key t cell features
features <- c("CD4", "IL7R", "GZMA","CCL4","IL17RB",
              "CD8A", "TCF7","GZMB", "GZMK", "IL13")
p <- prettyFeats(seu.obj = seu.obj, nrow = 2, ncol = 5, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 16, noLegend = T)
ggsave(paste("./output/", outName, "/",outName, "_fig4b.png", sep = ""), width = 12.5, height = 5)


### Fig supp 6a - feature plot of more key t cell features
features <- c("CXCR4", "SELL", "CCR7",
              "ITGB7", "ITGAE")
p <- prettyFeats(seu.obj = seu.obj, nrow = 2, ncol = 3, features = features, legInLine = T,
                 color = "black", order = F, pt.size = 0.0000001, title.size = 16, noLegend = F)
ggsave(paste("./output/", outName, "/",outName, "_supp6a.png", sep = ""), width = 9.5, height = 6)

### Fig extra - feature plot of more key t cell features
features <- c("CD3D", "CD3E", "CD3G","CD5")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 4, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 16, noLegend = F)
ggsave(paste("./output/", outName, "/",outName, "_extraFeats.png", sep = ""), width = 12.5, height = 3)


### Fig 4c - DEG between resident and Non-resident -- this generated supplemental_data_5
idents.1 = "IL7R_high"
idents.2 = "GZMA_high"
groupBy = "majorID_sub"
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = groupBy, idents.1 = idents.1, idents.2 = idents.2, bioRep = "name2",padj_cutoff = 0.05, lfcCut = 0.58, 
                      minCells = 25, outDir = paste0("./output/", outName, "/"), 
                      title = paste0(idents.1,"_vs_",idents.2), idents.1_NAME = idents.1, idents.2_NAME = idents.2, strict_lfc = T, 
                      returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 24, dwnCol = "#75149D", stblCol = "grey",upCol = "#FF007F", labSize = 3.5
                    )
p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in IL7R_high", leftLab = paste0("Up in GZMA_high"), arrowz = T,  rightCol = "#FF007F", leftCol = "#75149D") + labs(x = "log2(Fold change)") + NoLegend() + theme(panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 2),
                                      axis.line = element_blank())
ggsave(gsub(" ", "_",paste("./output/", outName, "/", outName, "_", idents.1, "_V_", idents.2, ".png", sep = "")), width = 7, height = 7)


#remane idents to split the Memory subtype into their cooresponding resdent/non-resident classification
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "GZMA_high", "1" = "GZMA_high", 
                                   "2" = "IL7R_high", "3" = "GZMA_high",
                                   "4" = "GZMA_high", "5" = "IL7R_high",
                                   "6" = "IL7R_high", "7" = "IL7R_high",
                                   "8" = "IL7R_high", "9" = "GZMA_high",
                                   "10" = "GZMA_high", "11" = "IL7R_high",
                                   "12" = "GZMA_high", "13" = "ILC2",
                                   "14" = "IL7R_high","15" = "IL7R_high")
                       )
seu.obj$majorID_sub_SplitMem <- Idents(seu.obj)

#remove the ILCs
seu.obj.sub <- subset(seu.obj, invert = T, subset = majorID_sub_SplitMem == "ILC2")
seu.obj.sub$majorID_sub_SplitMem <- droplevels(seu.obj.sub$majorID_sub_SplitMem)


### Fig extra - Evlauate cell frequency by majorID_sub_SplitMem -- results = trendy, but n.s.
freqy <- freqPlots(seu.obj.sub, method = 1, nrow= 1, groupBy = "majorID_sub_SplitMem", legTitle = "Cell source",refVal = "name2", showPval = F,
                   namez = unique(seu.obj$name2), 
                   colz = unique(seu.obj$colz)
                  )
freqy + ggpubr::stat_compare_means(method = "t.test",
                                   method.args = list(var.equal = F),
                                   aes(label = paste0("p = ", ..p.format..)), label.x.npc = "left", label.y.npc = 1,vjust = -1, size = 3)
ggsave(paste("./output/", outName, "/",outName, "_freqPlots_majorID_sub_SplitMem.png", sep = ""), width = 10, height = 5)


### Fig Supp 8c - look at GZMA_high/IL7R_high ratio
groupByList <- seu.obj.sub$name2
clusterList <- seu.obj.sub$majorID_sub_SplitMem

cluster_freq.table <- as.data.frame(table(groupByList, clusterList)) %>% melt()
cluster_freq.table <- cluster_freq.table[,-3]
colnames(cluster_freq.table) <- c("Sample", "ClusterID", "Count")
cluster_freq.table <- cluster_freq.table %>% dplyr::group_by(Sample) %>% mutate(pct = round(prop.table(Count),2))

ratio.df <- cluster_freq.table[ ,c(1:3)] %>% pivot_wider(names_from = "ClusterID", values_from = "Count") %>% mutate(trm_tinf_ratio = GZMA_high / IL7R_high)

ratio.df <- ratio.df[ ,c(1,4)] %>% mutate(cellSource = ifelse(grepl("CIE", Sample), "CIE", "Healthy"))
t.test(trm_tinf_ratio ~ cellSource, ratio.df)
ratio.df$cellSource <- factor(ratio.df$cellSource, levels = c("Healthy","CIE"))

car::leveneTest(trm_tinf_ratio ~ cellSource, data = ratio.df)

p <- ggplot(ratio.df, aes(y = trm_tinf_ratio, x = cellSource)) + 
    labs(x = NULL, y = "Ratio (GZMA_high/IL7R_high)") + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), 
          strip.background = element_rect(fill = NA, color = NA), 
          strip.text = element_text(face = "bold"), 
          axis.ticks.x = element_blank(), 
          axis.text = element_text(color = "black")         ) + 
    guides(fill = "none") + 
    geom_boxplot(aes(x = cellSource), alpha = 0.25, outlier.color = NA) + 
    geom_point(size = 2, position = position_jitter(width = 0.25),
               aes(x = cellSource, y = trm_tinf_ratio, color = Sample)) +
    labs(color = "Cell Source") + 
    ggpubr::stat_compare_means(method = "t.test", method.args = list(var.equal = T), aes(label = paste0("p = ", ..p.format..)), label.x.npc = "left", label.y.npc = 1, vjust = -1, size = 3) + 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    theme(panel.grid.major = element_line(color = "grey", size = 0.25),
          #legend.position = "none",
          text = element_text(size = 12) 
          ) + scale_color_manual(labels = unique(seu.obj$name2), values = unique(seu.obj$colz)) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_supp8c.png", sep = ""), width = 2.5, height = 3.5)


### Fig supp 6a: gsea of the DGE results
p <- plotGSEA(pwdTOgeneList = "./output/tcell/IL7R_high_vs_GZMA_high_all_genes.csv", category = "C5", 
              upCol = "#75149D", dwnCol = "#FF007F", size = 3.5) 

minVal <- -12
maxVal <- 12
pi <- p + scale_x_continuous(limits = c(minVal, maxVal), name = "Signed log10(padj)") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16)
    ) + 
    labs(
        title = "Gene ontology of T cell subsets",
        subtitle = "(IL7R_high vs GZMA_high)"
    ) + 
    annotate("segment", x = -0.1, 
             y = nrow(p$data)+1, 
             xend = minVal, 
             yend = nrow(p$data)+1, 
             lineend = "round", linejoin = "bevel", linetype ="solid", colour = "#75149D",
             size = 1, arrow = arrow(length = unit(0.1, "inches"))
            ) + 
    annotate(geom = "text", x = (minVal-0.1*1.5)/2-0.1*1.5, 
             y = nrow(p$data)+2,
             label = "Up in GZMA_high",
             hjust = 0.5,
             vjust = 1.5,
             size = 5) +
    annotate("segment", x = 0.1, 
             y = nrow(p$data)+1, 
             xend = maxVal,
             yend = nrow(p$data)+1,
             lineend = "round", linejoin = "bevel", linetype ="solid", colour = "#FF007F",
             size = 1, arrow = arrow(length = unit(0.1, "inches"))
            ) + 
    annotate(geom = "text", x = (maxVal-0.1*1.5)/2+0.1*1.5, 
             y = nrow(p$data)+2,
             label = "Up in IL7R_high",
             hjust = 0.5,
             vjust = 1.5,
             size = 5)
ggsave(paste("./output/", outName, "/", outName, "_supp6a.png", sep = ""), width = 7.6, height = 7)


### Fig extra - plot the top 20 DEGs in heatmap form
geneLists <- read.csv("./output/tcell/Non-resident_vs_Tissue_resident_all_genes.csv")

#extract top 20 feats for each direction of the conntrast
geneListUp <- geneLists %>% arrange(padj) %>% filter(log2FoldChange > 0) %>% .$gene
geneListDwn <- geneLists %>% arrange(padj) %>% filter(log2FoldChange < 0) %>% .$gene
feats_forHeat <- c(head(geneListUp,20), head(geneListDwn,20))

#subset the data
seu.obj.sub <- subset(seu.obj, majorID_sub_SplitMem == "Non-resident" | majorID_sub_SplitMem == "GZMA_high")
seu.obj.sub$type <- paste0(seu.obj.sub$majorID_sub_SplitMem,"--",seu.obj.sub$name2)

#extract metadata and data
metadata <- seu.obj.sub@meta.data
expression <- as.data.frame(t(seu.obj.sub@assays$RNA@data))
expression$anno_merge <- seu.obj.sub@meta.data[rownames(expression),]$type

#get cell type expression averages - do clus avg expression by sample
clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL     

#scale data
mat <- scale(as.matrix(clusAvg_expression))
mat <- mat[,feats_forHeat]

#prep annotations
samp <- unique(seu.obj$colz)
names(samp) <- unique(seu.obj$name2)
ha = HeatmapAnnotation(
    Sample = unlist(lapply(rownames(mat), function(x){strsplit(x,"--")[[1]][2]})),
    border = TRUE,
    col = list(Sample = samp)
)

#plot the data
png(file = paste0("./output/", outName, "/", outName, "_tcellHeat_degs.png"), width=3000, height=4000, res=400)
par(mfcol=c(1,1))         
ht <- Heatmap(t(mat), #name = "mat", #col = col_fun,
              name = "Scaled expression",
              cluster_rows = F,
              show_row_names=T,
              show_column_names=F,
              top_annotation = ha,
              col=viridis(option = "magma",100),
              cluster_columns = T
             )
draw(ht, padding = unit(c(2, 12, 2, 5), "mm"))
dev.off()


##################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin T cell analysis w/out ILC2s  ######## <<<<<<<<<<<<<<
##################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#remove ILC2s and recluster
seu.obj.sub <- subset(seu.obj, invert = T,
                  subset = 
                  clusterID_sub ==  "13")

table(seu.obj.sub$clusterID_sub)
table(seu.obj.sub$clusterID_2_1)
table(seu.obj.sub$orig.ident)

rm(seu.obj)
gc()

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "./output/s2/", subName = "240112_tcell_noILC_duod_h3c4_NoIntrons_2500", 
                     preSub = T, nfeatures = 2500, vars.to.regress = "percent.mt", saveRDS = F
                    )

#clustree to determine clus resolution
# seu.obj <- readRDS(file = "./output/s2/230913_tcell_duod_h3c4_NoIntrons_2500_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "240112_tcell_noILC_duod_h3c4_NoIntrons_2500", test_dims = c("40","35", "30"), algorithm = 3, prefix = "integrated_snn_res.")

#visulize the data & evaluate
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "240112_tcell_noILC_duod_h3c4_NoIntrons_2500", final.dims = 35, final.res = 0.6, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

### load in the obj
seu.obj <- readRDS("./output/s3/240112_tcell_noILC_duod_h3c4_NoIntrons_2500_res0.6_dims30_dist0.1_neigh10_S3.rds")
outName <- "tcell"
colz <- c("#A41DDC", "#C47AEA", "#DAACF2", "#75149D", "#FF007F", "#F77FBE", "#FF1E6D", "#9B1664", "#FF007F", "#F77FBE")
labColz <- c("white", "black", "black", "white", "black", "black", "black", "white", "black", "black")

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
seu.obj <- RenameIdents(seu.obj, c("0" = "CD8_eff (c0)", "1" = "CD8_eff (c0)", 
                                   "2" = "Tnaive (c1)", "3" = "CD8_TRM (c2)",
                                   "4" = "gdT_1 (c3)", "5" = "Tnaive (c1)",
                                   "6" = "NK_T (c4)", "7" = "Tnaive (c1)",
                                   "8" = "CD8_mem (c5)", "9" = "gdT_2 (c6)",
                                   "10.1" = "Tnaive (c1)", "10.2" = "CD8_eff (c0)",
                                   "11" = "CD8_eff (c0)", "12" = "Treg (c7)",
                                   "13" = "NK (c8)","14" = "T_IFN (c9)")
                       )
seu.obj$majorID_sub_noILC <- Idents(seu.obj)
seu.obj$majorID_sub_noILC <- factor(seu.obj$majorID_sub_noILC, levels = levels(seu.obj$majorID_sub_noILC)[c(4,7,1,3,6,2,8,10,5,9)])

Idents(seu.obj) <- "majorID_sub_big"
seu.obj <- RenameIdents(seu.obj, c("Non-resident" = "IL7R_high",
                                   "Tissue resident" = "GZMA_high",
                                   "Memory" = "CCL4_high"
                                  )
                       )
seu.obj$majorID_sub_big <- Idents(seu.obj)

#stash the numerical ID
clusterID_final <- table(seu.obj$majorID_sub_noILC) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
mutate(clusterID_final=row_number()-1) %>% arrange(clusterID_final) 

newID <- clusterID_final$clusterID_final
names(newID) <- clusterID_final$Var1
Idents(seu.obj) <- "majorID_sub_noILC"
seu.obj <- RenameIdents(seu.obj, newID)
table(Idents(seu.obj))
seu.obj$clusterID_final_noILC <- Idents(seu.obj)
table(seu.obj$majorID_sub_noILC, seu.obj$clusterID_final_noILC)
seu.obj$clusterID_final_noILC <- factor(seu.obj$clusterID_final_noILC, levels = levels(seu.obj$clusterID_final_noILC)[c(4,7,1,3,6,2,8,10,5,9)])


### supp data 7 - cell type gene signatures
#generate violin plots for each cluster
saveName <- "240128_tcell_noILC2_duod_h3c4_NoIntrons_2500"
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_sub_noILC", numOfFeats = 24, outName = saveName,
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )

#export surface marker data with FindAllMarkers
surface.markers <- read.csv("./surface_master.csv")[ ,c("UniProt.gene", "UniProt.description", "Surfaceome.Label", "Surfaceome.Label.Source")] %>% filter(!duplicated(UniProt.gene))
cluster.markers <- read.csv("./output/viln/tcell/240128_tcell_noILC2_duod_h3c4_NoIntrons_2500_gene_list.csv")
write.csv(cluster.markers[ ,c(7,8,2:6)] %>% left_join(surface.markers, by = c("gene" = "UniProt.gene")),
          file = "./output/supplementalData/supplemental_data_7.csv", row.names = F)


seu.obj$celltype.l1 <- seu.obj$majorID_sub_big
seu.obj$celltype.l2 <- seu.obj$majorID_sub_noILC
seu.obj$sample_name <- seu.obj$name2
### supp data - export data for cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "./output/cb_input/", 
               markers = "./output/supplementalData/supplemental_data_7.csv", 
               metaSlots = c("cluster","gene","avg_log2FC","p_val_adj", "UniProt.description", "Surfaceome.Label", "Surfaceome.Label.Source"),
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase",
                                                   "sample_name", "cellSource", "clusterID_sub", "celltype.l1", "celltype.l2"), 
               skipEXPR = T, test = F,
               feats = c("CD4", "CD8A", "GZMA", "GZMB", 
                         "CD40LG", "TCF7", "FASLG",
                         "THY1", "BIN1", "TNFRSF6B", "CTLA4",
                         "IL13", "IL17RB", "NCR3","F2RL3")
                          )

seu.obj <- cleanMeta(seu.obj = seu.obj, 
                     metaSlot_keep = c(
                         "orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", "majorID", 
                         "nCount_SCT", "nFeature_SCT", "clusterID",
                         "colz", "sample_name", "cellSource", "clusterID_sub", "celltype.l1", "celltype.l2"
                     )
                    )
saveRDS(seu.obj, file = "./output/s3/Tcell_duod_annotated.rds")



### Fig supp 8a - feature plot of key t cell features
features <- c("CD4", "IL7R", "GZMA","CCL4",
              "CD8A", "TCF7","GZMB", "GZMK", 
              "CTLA4", "NCR3", "OAS1", "IFNG"
             )
p <- prettyFeats(seu.obj = seu.obj, nrow = 3, ncol = 4, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 16, noLegend = F)
ggsave(paste("./output/", outName, "/",outName, "_supp8a.png", sep = ""), width = 12.5, height = 9)


### Fig 4d - tcell colorized by annotation clustering
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_final_noILC",
              pt.size = 0.25,
              cols = colz,
              label = T,
              label.box = T
             ) + NoLegend()
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.8, labCol = labColz, smallAxes = T)
ggsave(paste0("./output/", outName, "/", outName, "_fig4d.png"), width = 7, height = 7)


### Fig extra - tcell colorized by annotation clustering
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "majorID_sub_noILC",
              pt.size = 0.25,
              cols = colz,
              label = T,
              label.box = T
             )
p <- formatUMAP(plot = pi) + theme(axis.title = element_blank(),
                             panel.border = element_blank(),
                             plot.margin = unit(c(-7, -7, -7, -7), "pt")
                            )
ggsave(paste0("./output/", outName, "/", outName, "_lab_umap.png"), width = 9, height = 7)


### Fig supp 8a - sankey plot to show how clustering changed
#use custom function then modify the plot (TO DO: fix the custom fuction versitility)
p <- sankeyPlot(seu_obj = seu.obj, new.ident = "clusterID_final_noILC", old.ident = "majorID_sub_big", old.colorz = gg_color_hue(3),
                       new.colorz =  gg_color_hue(13), old.labCol = "black", new.labCol = "black", flowCol = "grey"
                    )

#extract and reorder nodes
data <- p$data %>% as.data.frame()
data$node <- factor(data$node,levels = c("GZMA_high", "CCL4_high", "IL7R_high", "c0", "c3", "c6", "c2", "c5", "c1", "c7", "c4", "c8", "c9"))
data$next_node <- factor(data$next_node,levels = c("c0", "c3", "c6", "c2", "c5", "c1", "c7", "c4", "c8", "c9"))

#color order is alphabetical for some reason... replot with correct colorization
p <- ggplot(data, aes(x = x,
                      next_x = next_x,
                      node = node,
                      next_node = next_node,
                      fill = factor(node),
                      label = node)) + 
geom_sankey() +
geom_sankey(flow.alpha = 0.5, node.color = 1) + 
geom_sankey_label(size = 3, color = "black", fill = "gray85") +
theme_sankey(base_size = 16) +
scale_fill_manual(values = c("#FF1E6D", "#C47AEA", "#EF98CB", colz[c(3,6,4,1,9,5,2,7,10,8)])) +
theme(legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(t = 7, r = 7, b = 7, l = 28, unit = "pt")) + scale_x_discrete(expand = c(0, 0)) + coord_cartesian(clip = "off")
ggsave(paste0("./output/", outName, "/", outName, "_supp8a.png"), width = 3, height = 7)


### Fig extra - dot plot of key t cell features
features <- c("ENSCAFG00000041010","ENSCAFG00000011098","ENSCAFG00000030206",
              "GZMA","GZMB",
              "CCL4","GZMK",
              "IL7R","VIM",
              "CTLA4","FOXP3","TNFRSF4","TNFRSF18",
              "IFI44","OAS1","ISG15",
              "TCF7","CRTAM","CXCR4","NCR3","F2RL3","KLRB1","IL12RB2")
pi <- majorDot(seu.obj = seu.obj, groupBy = "majorID_sub_noILC",
               features = features
              ) + theme(axis.title = element_blank(),
                        axis.text = element_text(size = 12))
ggsave(paste("./output/", outName, "/", outName, "_noILC_MajorDots.png", sep = ""), width = 9, height = 5)


### Fig supp 8b - transfer TRDC data from ROS data
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

p <- prettyFeats(seu.obj = seu.obj.Hsub, nrow = 2, ncol = 1, features = features, color = "black", order = F) 
ggsave(paste("./output/", outName, "/", outName, "_supp8b.png", sep = ""), width = 5, height = 8)


### Fig supp 5e - umap by sample
Idents(seu.obj) <- "cellSource"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data$cellSource)))
pi <- DimPlot(seu.obj.ds, 
              reduction = "umap", 
              group.by = "name2",
              cols = unique(seu.obj.ds$colz), #check colorization is correct
              pt.size = 0.25,
              label = FALSE,
              shuffle = TRUE
)
pi <- pi + labs(colour="") + theme(axis.title = element_blank(),
                                               panel.border = element_blank(),
                                                  legend.position = "top", 
                                              legend.justification = "center",
                                              legend.direction = "horizontal",legend.title=element_text(size=12)) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))
p <- formatUMAP(pi, smallAxes = T)
ggsave(paste("./output/", outName, "/", outName, "_sup5e.png", sep = ""), width =7, height = 7)


### Fig supp 5f - evlauate cell frequency by majorID_sub_noILC
freqy <- freqPlots(seu.obj, method = 1, nrow = 2, groupBy = "majorID_sub_noILC", legTitle = "Cell source",refVal = "name2", showPval = T,
                   namez = unique(seu.obj$name2), 
                   colz = unique(seu.obj$colz)
                  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + NoLegend()

ggsave(paste("./output/", outName, "/",outName, "_sup5f.png", sep = ""), width = 12, height = 6)

#confirm variability is consistent between CIE and H
res.ftest <- lapply(levels(freqy$data$majorID_sub), function (x){
    data.df <- freqy$data[freqy$data$majorID_sub ==  x, ]
    car::leveneTest(freq ~ cellSource, data = data.df)
})


### Fig 4e - heatmap by cluster for each sample with CD4/CD8 proportion
seu.obj$majorID_sub_noILC <- factor(seu.obj$majorID_sub_noILC, levels = levels(seu.obj$majorID_sub_noILC)[c(1:5,8,6,7,9,10)])
colz <- c("#A41DDC", "#C47AEA", "#DAACF2", "#75149D", "#FF007F", "#F77FBE", "#FF1E6D", "#9B1664", "#FF007F", "#F77FBE")[c(1:5,8,6,7,9,10)]
seu.obj$type <- paste0(seu.obj$majorID_sub_noILC,"-",seu.obj$name2)
keep <- names(table(seu.obj$type))[table(seu.obj$type) >= 5]
Idents(seu.obj) <- "type"
seu.obj <- subset(seu.obj, idents = keep)

#calc the proportion of cells cd4 or cd8 positive
exp.df <- FetchData(object = seu.obj, vars = c('cellSource', 'type','name2', 'majorID_sub_noILC', 'CD4','CD8A'))
exp.df <- exp.df %>% group_by(name2,majorID_sub_noILC) %>% mutate(cd4_pos = ifelse(CD4 > 0 & CD8A == 0, 1, 0),
                                     cd8a_pos = ifelse(CD8A > 0 & CD4 == 0, 1, 0),
                                     total_cellz = n()
                                    ) %>% summarize(CD4 = sum(cd4_pos)/total_cellz,
                                                    CD8A = sum(cd8a_pos)/total_cellz) %>% distinct() %>% mutate(paste0(majorID_sub_noILC,"-",name2)) %>% as.data.frame()
exp.df$majorID_sub_noILC <- factor(exp.df$majorID_sub_noILC, levels = levels(seu.obj$majorID_sub_noILC))
exp.df$odor <- as.numeric(exp.df$majorID_sub_noILC)
exp.df <- exp.df %>% arrange(odor)

#extract metadata and data
metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$RNA@data)) #use log noralized count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type

#get cell type expression averages - do clus avg expression by sample
clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL

df <- rownames(clusAvg_expression) %>% as.data.frame()
colnames(df) <- "type"
df$clus <- unlist(lapply(df$type, function(x){strsplit(x,"-")[[1]][1]}))
df$clus <- factor(df$clus, levels = levels(seu.obj$majorID_sub_noILC))
df$odor <- as.numeric(df$clus)
df <- df %>% arrange(odor)

#load in feats that define
sig.df <- read.csv("./output/viln/tcell/240118_tcell_noILC_h3c4_NoIntrons_2500_gene_list.csv")
sig.df <- sig.df %>% filter(p_val_adj < 0.01)
sig.df$cluster <- factor(sig.df$cluster, levels = levels(seu.obj$majorID_sub_noILC))
sig.df$odor <- as.numeric(sig.df$cluster)

#extract labels to plot
lab.df <- sig.df[!grepl("ENSCAF", sig.df$gene), ]
text_list <- rev(split(lab.df$gene, lab.df$cluster))
text_list <- lapply(text_list, function(x){c(paste0(paste(x[1:3], collapse = ", "),","), paste(x[4:6], collapse = ", "))})

#bring over the one ENSCAF needed
text_list$`gdT_2 (c6)` <- sig.df[sig.df$cluster == "gdT_2 (c6)", ]$gene

#finish ordering
sig.df <- sig.df[!duplicated(sig.df[,"gene"]), ]
sig.df <- sig.df %>% arrange(odor)
geneOrder <- rev(sig.df$gene)

#filter matrix for feats that define and scale by row
clusAvg_expression <- clusAvg_expression[ ,colnames(clusAvg_expression) %in% sig.df$gene]
mat_scaled <- t(apply(t(clusAvg_expression), 1, scale))
colnames(mat_scaled) <- rownames(clusAvg_expression)
length(geneOrder) == nrow(mat_scaled)
mat_scaled <- mat_scaled[ ,match(df$type, colnames(mat_scaled))]
mat_scaled <- mat_scaled[match(geneOrder, rownames(mat_scaled)), ]        

#set annotations
samp <- unique(seu.obj$colz)
names(samp) <- unique(seu.obj$name2)
clus <- colz
names(clus) <- levels(seu.obj$majorID_sub_noILC)
CD4_col <- colorRamp2(c(0, max(exp.df$CD4)), c("white", "#FF007F"))
CD8_col <- colorRamp2(c(0, max(exp.df$CD8)), c("white", "purple"))
heat_col <- viridis(option = "magma",100)
ha <- HeatmapAnnotation(
    Sample = unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"-")[[1]][2]})),
    Cluster = unlist(lapply(colnames(mat_scaled), function(x){strsplit(x,"-")[[1]][1]})),
    CD4 = exp.df$CD4,
    CD8 = exp.df$CD8,
    border = TRUE,
    col = list(Sample = samp, Cluster = clus, CD4 = CD4_col, CD8 = CD8_col),
    show_legend = c(TRUE, FALSE, FALSE, FALSE),
    annotation_legend_param = list(
        Sample = list(
            nrow = 1
        )
    )
)

#make lgnd seperately
lgd1 <- Legend(labels = sig.df[!duplicated(sig.df[,"cluster"]), ]$cluster, legend_gp = gpar(fill = clus), title = "Clusters", 
    nrow = 2, gap = unit(2.5, "cm"))
lgd2 <- Legend(col_fun = CD4_col,
              at = c(0, round((max(exp.df$CD4)/2),2), round(max(exp.df$CD4),2)),
              title = "Proportion CD4+",
              direction = "horizontal")
lgd3 <- Legend(col_fun = CD8_col,
              at = c(0, round((max(exp.df$CD8)/2),2), round(max(exp.df$CD8),2)),
              title = "Proportion CD8+",
              direction = "horizontal")
ra <- rowAnnotation(foo = anno_empty(border = FALSE, width = max_text_width(unlist(text_list)) + unit(4, "mm")))
pd <- packLegend(lgd1, lgd2, lgd3, max_width = unit(10, "cm"), 
    direction = "horizontal", column_gap = unit(10, "mm"), row_gap = unit(0.5, "cm"))

#plot the data
ht <- Heatmap(
    mat_scaled,
    name = "mat",
    cluster_rows = F,
    row_title_gp = gpar(fontsize = 24),
    show_row_names=F,
    col=heat_col,
    cluster_columns = F,
    top_annotation = ha,
    right_annotation = ra,
    show_column_names = F,
    column_split = df$clus,
    row_split = factor(rev(sig.df$cluster), levels = rev(sig.df[!duplicated(sig.df[,"cluster"]), ]$cluster)),
    cluster_row_slices=F,
    row_title = NULL,
    column_title = NULL,
    heatmap_legend_param = list(
            title = "Scaled expression",
            direction = "horizontal"
        )
)

#save the plot
png(file = paste0("./output/", outName, "/", outName, "_fig3e.png"), width=3500, height=4000, res=400)
par(mfcol=c(1,1))    
draw(ht, padding = unit(c(2, 2, 2, 2), "mm"), annotation_legend_list = pd, merge_legend = TRUE, heatmap_legend_side = "bottom", 
    annotation_legend_side = "bottom")

for(i in 1:length(rev(sig.df[!duplicated(sig.df[,"cluster"]), ]$cluster))) {
    decorate_annotation("foo", slice = i, {
        grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(1, "mm"), just = "left",
                  gp = gpar(fontsize = 10))
    })
}

for(i in 1:length(rev(sig.df[!duplicated(sig.df[,"cluster"]), ]$cluster))) {
    decorate_annotation("Cluster", slice = i, {
        grid.text(unlist(lapply(levels(df$clus), function(x){gsub("\\)|\\(", "",strsplit(x," ")[[1]][2])}))[i], just = "center", gp = gpar(col = labColz[c(1:5,8,6,7,9,10)][i]))
    })
}
dev.off()


################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end T cell analysis w/out ILC2s  ######## <<<<<<<<<<<<<<
################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<