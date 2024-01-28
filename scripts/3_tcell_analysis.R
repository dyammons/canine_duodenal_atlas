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


#remove susspected low quality cluster
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
seu.obj <- readRDS("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/s3/230913_tcell_duod_h3c4_NoIntrons_2500_res0.6_dims35_dist0.3_neigh30_S3.rds")
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


### Fig extra - check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/",outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


### Fig extra - run singleR
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
ggsave(paste0("./output/", outName, "/", outName, "_supp2c.png"), width = 7, height = 7)


### supp data - cell type gene signatures
#generate violin plots for each cluster
saveName <- "231022_tcell_duod_h3c4_NoIntrons_2500"
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_sub", numOfFeats = 24, outName = saveName,
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )

#export surface marker data with FindAllMarkers
surface.markers <- read.csv("./surface_master.csv")[ ,c("UniProt.gene", "UniProt.description", "Surfaceome.Label", "Surfaceome.Label.Source")] %>% filter(!duplicated(UniProt.gene))
cluster.markers <- read.csv("./output/viln/tcell/231022_tcell_duod_h3c4_NoIntrons_2500_gene_list.csv")
write.csv(cluster.markers[ ,c(7,8,2:6)] %>% left_join(surface.markers, by = c("gene" = "UniProt.gene")),
          file = "./output/supplementalData/tcell_markers.csv", row.names = F)


### supp data - export data for cell browser
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


### Fig 2a - tcell colorized by annotation clustering
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_final",
              pt.size = 0.25,
              cols = colz,
              label = F,
              label.box = F
             )
p <- formatUMAP(plot = pi) + theme(axis.title = element_blank(),
                             panel.border = element_blank(),
                             plot.margin = unit(c(-7, -7, -7, -7), "pt")
                            ) + NoLegend()

ggsave(paste0("./output/", outName, "/", outName, "_fig2a.png"), width = 7, height = 7)


### Fig 2b - dot plot of key t cell features
features <- c("CD4", "IL7R", "GZMA","CCL4","IL17RB",
              "CD8A", "TCF7","GZMB", "GZMK", "IL13")

### Fig extra - plot on featPlots
p <- prettyFeats(seu.obj = seu.obj, nrow = 2, ncol = 5, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 16, noLegend = T)
ggsave(paste("./output/", outName, "/",outName, "_tcell_featPlots.png", sep = ""), width = 12.5, height = 5)

#remane idents to split the Memory subtype into their cooresponding resdent/non-resident classification
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "Tissue resident", "1" = "Tissue resident", 
                                   "2" = "Non-resident", "3" = "Tissue resident",
                                   "4" = "Tissue resident", "5" = "Non-resident",
                                   "6" = "Non-resident", "7" = "Non-resident",
                                   "8" = "Non-resident", "9" = "Tissue resident",
                                   "10" = "Tissue resident", "11" = "Non-resident",
                                   "12" = "Tissue resident", "13" = "ILC2",
                                   "14" = "Non-resident","15" = "Non-resident")
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
ggsave(paste("./output/", outName, "/",outName, "_freqPlots.png", sep = ""), width = 10, height = 5)


### Look at Trm/Tnaive ratio
groupByList <- seu.obj.sub$name2
clusterList <- seu.obj.sub$majorID_sub_SplitMem

cluster_freq.table <- as.data.frame(table(groupByList, clusterList)) %>% melt()
cluster_freq.table <- cluster_freq.table[,-3]
colnames(cluster_freq.table) <- c("Sample", "ClusterID", "Count")
cluster_freq.table <- cluster_freq.table %>% dplyr::group_by(Sample) %>% mutate(pct = round(prop.table(Count),2))

ratio.df <- cluster_freq.table[ ,c(1:3)] %>% pivot_wider(names_from = "ClusterID", values_from = "Count") %>% mutate(trm_tinf_ratio = `Tissue resident`/`Non-resident`)

ratio.df <- ratio.df[ ,c(1,4)] %>% mutate(cellSource = ifelse(grepl("CIE", Sample), "CIE", "Healthy"))
t.test(trm_tinf_ratio ~ cellSource, ratio.df)
ratio.df$cellSource <- factor(ratio.df$cellSource, levels = c("Healthy","CIE"))

car::leveneTest(trm_tinf_ratio ~ cellSource, data = ratio.df)

p <- ggplot(ratio.df, aes(y = trm_tinf_ratio, x = cellSource)) + 
    labs(x = NULL, y = "Ratio (Tissue resident/Non-resident)") + 
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
ggsave(paste("./output/", outName, "/", outName, "_RatioBoxplot.png", sep = ""), width = 2.5, height = 3.5)


### Fig 1d - stacked bar graph by major cell types
Idents(seu.obj) <- "cellSource"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data$cellSource)))

p <- stackedBar(seu.obj = seu.obj.ds, downSampleBy = "cellSource", groupBy = "name2", clusters = "majorID_sub_SplitMem") +
scale_fill_manual(labels = unique(seu.obj$name2), 
                  values = unique(seu.obj$colz)) + theme(axis.title.y = element_blank(),
                                                         axis.title.x = element_text(size = 14),
                                                         axis.text = element_text(size = 12)
                                                        ) 
ggsave(paste("./output/", outName,"/", outName, "_fig1d.png", sep = ""), width =7, height = 4)


### Fig 2c - DEG between resident and Non-resident
idents.1 = "Non-resident"
idents.2 = "Tissue resident"
groupBy = "majorID_sub_SplitMem"
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = groupBy, idents.1 = idents.1, idents.2 = idents.2, bioRep = "name2",padj_cutoff = 0.05, lfcCut = 0.58, 
                      minCells = 25, outDir = paste0("./output/", outName, "/"), 
                      title = paste0(idents.1,"_vs_",idents.2), idents.1_NAME = idents.1, idents.2_NAME = idents.2, 
                      returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 24, dwnCol = "#75149D", stblCol = "grey",upCol = "#FF007F", labSize = 3.5
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = paste0("Up in ",idents.1), leftLab = paste0("Up in ",idents.2), arrowz = T,  rightCol = "#FF007F", leftCol = "#75149D") + labs(x = paste0("log2(FC)", idents.1, " vs ",idents.2)) + NoLegend() + theme(panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 2),
                                      axis.line = element_blank())
ggsave(gsub(" ", "_",paste("./output/", outName, "/", outName, "_", idents.1, "_V_", idents.2, ".png", sep = "")), width = 7, height = 7)

#colz <- c("#A41DDC", "#C47AEA", "#DAACF2", "#75149D", "#FF007F", "#F77FBE", "#FF1E6D", "#9B1664", "#FF007F", "#F77FBE")

### Fig supp: gsea of the DGE results
p <- plotGSEA(pwdTOgeneList = "./output/tcell/Non-resident_vs_Tissue_resident_all_genes.csv", category = "C5", 
              upCol = "blue", dwnCol = "red", size = 4) 

minVal <- -10
maxVal <- 10

pi <- p + scale_x_continuous(limits = c(minVal, maxVal), name = "Signed log10(padj)") + 
    theme(axis.title=element_text(size = 16)) + 

    annotate("segment", x = -0.1, 
             y = nrow(p$data)+1, 
             xend = minVal, 
             yend = nrow(p$data)+1, 
             lineend = "round", linejoin = "bevel", linetype ="solid", colour = "blue",
             size = 1, arrow = arrow(length = unit(0.1, "inches"))
            ) + 

    annotate(geom = "text", x = (minVal-0.1*1.5)/2-0.1*1.5, 
             y = nrow(p$data)+2,
             label = "Up in Non-resident",
             hjust = 0.5,
             vjust = 1.5,
             size = 5) +

    annotate("segment", x = 0.1, 
             y = nrow(p$data)+1, 
             xend = maxVal,
             yend = nrow(p$data)+1,
             lineend = "round", linejoin = "bevel", linetype ="solid", colour = "red",
             size = 1, arrow = arrow(length = unit(0.1, "inches"))
            ) + 

    annotate(geom = "text", x = (maxVal-0.1*1.5)/2+0.1*1.5, 
             y = nrow(p$data)+2,
             label = "Up in tissue resident",
             hjust = 0.5,
             vjust = 1.5,
             size = 5)


ggsave(paste("./output/", outName, "/", outName, "_fig2e.png", sep = ""), width = 7.6, height = 5)


### Fig extra - plot the top 20 DEGs in heatmap form
geneLists <- read.csv("./output/tcell/Non-resident_vs_Tissue_resident_all_genes.csv")

#extract top 20 feats for each direction of the conntrast
geneListUp <- geneLists %>% arrange(padj) %>% filter(log2FoldChange > 0) %>% .$gene
geneListDwn <- geneLists %>% arrange(padj) %>% filter(log2FoldChange < 0) %>% .$gene
feats_forHeat <- c(head(geneListUp,20), head(geneListDwn,20))

#subset the data
seu.obj.sub <- subset(seu.obj, majorID_sub_SplitMem == "Non-resident" | majorID_sub_SplitMem == "Tissue resident")
seu.obj.sub$type <- paste0(seu.obj.sub$majorID_sub_SplitMem,"--",seu.obj.sub$name2)

#extract metadata and data
metadata <- seu.obj.sub@meta.data
expression <- as.data.frame(t(seu.obj.sub@assays$RNA@data)) #use log noralized count
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
png(file = paste0("./output/", outName, "/", outName, "_tcellHeat.png"), width=3000, height=4000, res=400)
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


### Fig 2d - tcell colorized by annotation clustering
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_final_noILC",
              pt.size = 0.25,
              cols = colz,
              label = T,
              label.box = T
             )
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.8, labCol = labColz, smallAxes = T)
ggsave(paste0("./output/", outName, "/", outName, "_fig2d.png"), width = 7, height = 7)


### Fig supp 5b - tcell colorized by annotation clustering
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
ggsave(paste0("./output/", outName, "/", outName, "_fig2d.png"), width = 9, height = 7)



### Supp data: generate feats that define each cell type
saveName <- "240118_tcell_noILC_h3c4_NoIntrons_2500"
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_sub_noILC", numOfFeats = 24, outName = saveName,
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )


### Fig supp - sankey plot to show how clustering changed
#use custom function then modify the plot (TO DO: fix the custom fuction versitility)
p <- sankeyPlot(seu_obj = seu.obj, new.ident = "clusterID_final_noILC", old.ident = "majorID_sub_big", old.colorz = gg_color_hue(3),
                       new.colorz =  gg_color_hue(13), old.labCol = "black", new.labCol = "black", flowCol = "grey"
                    )

#extract and reorder nodes
data <- p$data %>% as.data.frame()
data$node <- factor(data$node,levels = c("Tissue resident", "Memory", "Non-resident", "c0", "c3", "c6", "c2", "c5", "c1", "c7", "c4", "c8", "c9"))
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
scale_fill_manual(values = c("#FF1E6D", "#EF98CB", "#C47AEA", colz[c(3,6,4,1,9,5,2,7,10,8)])) +
theme(legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(t = 7, r = 7, b = 7, l = 28, unit = "pt")) + scale_x_discrete(expand = c(0, 0)) + coord_cartesian(clip = "off")
ggsave(paste0("./output/", outName, "/", outName, "_fig2a.png"), width = 3, height = 7)


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



### Fig supp - evlauate cell frequency by majorID_sub_noILC
freqy <- freqPlots(seu.obj, method = 1, nrow = 2, groupBy = "majorID_sub_noILC", legTitle = "Cell source",refVal = "name2", showPval = T,
                   namez = unique(seu.obj$name2), 
                   colz = unique(seu.obj$colz)
                  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + NoLegend()

ggsave(paste("./output/", outName, "/",outName, "_fig2c.png", sep = ""), width = 12, height = 6)

#confirm variability is consistent between CIE and H
res.ftest <- lapply(levels(freqy$data$majorID_sub), function (x){
    data.df <- freqy$data[freqy$data$majorID_sub ==  x, ]
    car::leveneTest(freq ~ cellSource, data = data.df)
})


### Fig supp 2c - umap by sample
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

### Fig 2e - heatmap by cluster for each sample with CD4/CD8 proportion
library(circlize)
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
sig.df <- read.csv("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/viln/tcell/240118_tcell_noILC_h3c4_NoIntrons_2500_gene_list.csv")
sig.df <- sig.df %>% filter(p_val_adj < 0.01)
sig.df$cluster <- factor(sig.df$cluster, levels = levels(seu.obj$majorID_sub_noILC))
sig.df$odor <- as.numeric(sig.df$cluster)

#extract labels to plot
lab.df <- sig.df[!grepl("ENSCAF", sig.df$gene), ]
text_list <- rev(split(lab.df$gene, lab.df$cluster))
text_list <- lapply(text_list, function(x){c(paste0(paste(x[1:3], collapse = ", "),","), paste(x[4:6], collapse = ", "))})

#bring over the one ENSCAF needed
text_list$`gdT_2 (c6)` <- sig.df[sig.df$cluster == "gdT_2 (c6)", ]$gene

sig.df <- sig.df[!duplicated(sig.df[,"gene"]), ]
sig.df <- sig.df %>% arrange(odor)

#order the rows
geneOrder <- rev(sig.df$gene)

#filter matrix for feats that define and scale by row
clusAvg_expression <- clusAvg_expression[ ,colnames(clusAvg_expression) %in% sig.df$gene]
mat_scaled <- t(apply(t(clusAvg_expression), 1, scale))
# mat_scaled <- FastRowScale(t(as.matrix(clusAvg_expression)), center = TRUE, scale = TRUE, scale_max = 5)

colnames(mat_scaled) <- rownames(clusAvg_expression)

length(geneOrder) == nrow(mat_scaled)
mat_scaled <- mat_scaled[ ,match(df$type, colnames(mat_scaled))]
mat_scaled <- mat_scaled[match(geneOrder, rownames(mat_scaled)), ]        

# colnames(mat_scaled)
# [c(39:44,15:21,34:38,8:14,1:7)]

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
    mat_scaled, #name = "mat", #col = col_fun,
    name = "mat",
    cluster_rows = F,
#   row_title = "Genes",
    row_title_gp = gpar(fontsize = 24),
    show_row_names=F,
    col=heat_col,
    cluster_columns = F,
    top_annotation = ha,
    right_annotation = ra,
    show_column_names = F,
    column_split = df$clus,
    row_split = factor(rev(sig.df$cluster), levels = rev(sig.df[!duplicated(sig.df[,"cluster"]), ]$cluster)),
#     show_heatmap_legend = F,
#   row_km = length(rev(sig.df[!duplicated(sig.df[,"cluster"]), ]$cluster)),
    cluster_row_slices=F,
    row_title = NULL,
    column_title = NULL,
    heatmap_legend_param = list(
            title = "Scaled expression",
            direction = "horizontal"
        )
#   column_title = "Cell types",
#   column_title_gp = gpar(fontsize = 24),
#   column_title_side = "bottom",
#   column_names_rot = 45,
)

png(file = paste0("./output/", outName, "/", outName, "_tcellHeat.png"), width=3500, height=4000, res=400)
par(mfcol=c(1,1))    
# draw(ht, padding = unit(c(2, 12, 2, 5), "mm"), heatmap_legend_side = "top")
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








### Fig extra - plot by lower res clusters
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "majorID_sub_big",
              pt.size = 0.25,
              cols = colz.df$V1,
              label = F,
              label.box = F
             )
p <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste0("./output/", outName, "/", outName, "_major.png"), width = 7, height = 7)


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
                        axis.text = element_text(size = 12)) #+ geom_tile(aes(fill = c(1,1,1,2,2,2,2,2,2,3), x = 0), size = 1, show.legend = FALSE)
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
              shuffle = F
)
pi <- formatUMAP(plot = pi) + NoLegend() + theme(plot.title = element_text(size = 18, vjust = 1),
                                                 axis.title = element_blank(),
                                                 panel.border = element_blank()) + ggtitle("Canine PBMC reference mapping")
ggsave(paste("./output/", outName,"/",outName, "_umap_Predicted.png", sep = ""), width = 7, height = 7)


### Fig - reference map using gut atlas data
reference <- MuDataSeurat::ReadH5AD("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/Tcell_log_counts02_v2.h5ad")
#download reference from https://www.gutcellatlas.org/#datasets

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
                                                 panel.border = element_blank()) + ggtitle("Human T cell reference mapping")
ggsave(paste("./output/", outName,"/",outName, "_umap_Predicted_gutAtlas.png", sep = ""), width = 7, height = 7)

#wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_final_dataSet_H.rds.gz #just healthy
seu.pbmc <- readRDS("GSE225599_final_dataSet_H.rds")
Idents(seu.pbmc) <- "celltype.l1"
seu.pbmc <- subset(seu.pbmc, subset = celltype.l1 == "CD4 T cell" | celltype.l1 == "CD8/NK cell" | celltype.l1 == "gd T cell" | celltype.l1 == "DN T cell")
seu.pbmc.list <- SplitObject(seu.pbmc, split.by = "orig.ident")
seu.duod.list <- SplitObject(seu.obj, split.by = "orig.ident")

#integrate the two datasets
seu.obj <- indReClus(seu.list = c(seu.pbmc.list,seu.duod.list), outDir = "./output/s2/", subName = "pbmc_duod", preSub = T, nfeatures = 2000, saveRDS = F, returnObj = T, 
                      vars.to.regress = "percent.mt"
                       )

#visulize the data & evaluate
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "240111_tcell_pbmc_NoIntrons_2000", final.dims = 40, final.res = 0.6, stashID = "clusterID_sub_pbmc", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = F,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

seu.obj <- readRDS("./output/s3/240111_tcell_pbmc_NoIntrons_2000_res0.2_dims40_dist0.3_neigh30_S3.rds")

seu.obj$tissueSource <- ifelse(grepl("cfam", seu.obj$orig.ident), "pbmc", "duod") 
table(seu.obj$clusterID_sub_pbmc, seu.obj$tissueSource)

seu.obj$ct_merge <- as.factor(ifelse(grepl("cfam", seu.obj$orig.ident), as.character(seu.obj$celltype.l2), as.character(seu.obj$majorID_sub)))

Idents(seu.obj) <- "tissueSource"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data$tissueSource)))
pi <- DimPlot(seu.obj.ds, 
              reduction = "umap", 
              group.by = "ct_merge",
              #cols = levels(seu.obj.ds$colz), #check colorization is correct
              split.by = "tissueSource",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F
)
pi <- formatUMAP(plot = pi) + NoLegend() + theme(plot.title = element_text(size = 18, vjust = 1),
                                                 axis.title = element_blank(),
                                                 panel.border = element_blank()) + ggtitle("Canine PBMC + duod T cells")
ggsave(paste("./output/", outName,"/",outName, "_umap_int.png", sep = ""), width = 14, height = 7)


pi <- DimPlot(seu.obj.ds, 
              reduction = "umap", 
              group.by = "ct_merge",
              #cols = levels(seu.obj.ds$colz), #check colorization is correct
              split.by = "ct_merge",
              ncol = 5,
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F
)
pi <- formatUMAP(plot = pi) + NoLegend() + theme(plot.title = element_text(size = 18, vjust = 1),
                                                 axis.title = element_blank(),
                                                 panel.border = element_blank()) + ggtitle("Canine PBMC + duod T cells")
ggsave(paste("./output/", outName,"/",outName, "_umap_int.png", sep = ""), width = 18, height = 18)



### Fig extra: plot enrichment scores
ecLists <- read.csv("gut_ecTerms.csv", header = T)
ecLists <- ecLists[ecLists$lineage == "T_NK_cells", ]

modulez <- split(ecLists$genes, ecLists$cluster)
modulez <- modulez[unname(unlist(lapply(unlist(lapply(modulez, length)), function(x){ifelse(x >= 10, TRUE, FALSE)})))]
names(modulez) <- paste0(names(modulez),"_SIG")

seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")

names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)
features <- names(modulez)

p <- majorDot(seu.obj = seu.obj, groupBy = "majorID_sub",
                     features = features
                    ) + theme(legend.position = "bottom",
                              axis.title.y = element_blank(),
                              plot.margin = margin(7, 7, 0, 100, "pt")) + scale_y_discrete(position = "right") + guides(size = guide_legend(nrow = 2, byrow = F, title = 'Percent\nenriched')) + guides(color = guide_colorbar(title = 'Scaled\nenrichment score')) 
ggsave(paste("./output/", outName, "/", outName, "_gut_modScores.png", sep = ""), width = 8, height = 6)



### Fig extra: plot enrichment scores
# ecLists <- read.csv("./output/allCells/n3n4/linDEG/majorID_pertyName_T cellgeneList.csv", header = T)
# ecLists <- ecLists[ecLists$avg_log2FC > 0, ]

# modulez <- split(ecLists$X, ecLists$cellType)
# modulez <- modulez[unname(unlist(lapply(unlist(lapply(modulez, length)), function(x){ifelse(x >= 10, TRUE, FALSE)})))]
modulez <- list("CD8_EM_SIG" = c("CCL5", "GZMH", "CD8A", "TRAC", "KLRD1", "NKG7", "GZMK", "CST7", "CD8B", "TRGC2"))

seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")

names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)
features <- names(modulez)

p <- majorDot(seu.obj = seu.obj, groupBy = "majorID_sub",
                     features = features
                    ) + theme(legend.position = "bottom",
                              axis.title.y = element_blank(),
                              plot.margin = margin(7, 7, 0, 100, "pt")) + scale_y_discrete(position = "right") + guides(size = guide_legend(nrow = 2, byrow = F, title = 'Percent\nenriched')) + guides(color = guide_colorbar(title = 'Scaled\nenrichment score')) 
ggsave(paste("./output/", outName, "/", outName, "_gut_modScores.png", sep = ""), width = 8, height = 6)


features <- c("CD8_EM_SIG")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 1, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/",outName, "_QC_feats.png", sep = ""), width = 4, height = 4)





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

seu.obj.sub <- subset(seu.obj, invert = T, subset = majorID_sub_big == "ILC2")
seu.obj.sub$majorID_sub_big <- droplevels(seu.obj.sub$majorID_sub_big)
### Fig extra - Evlauate cell frequency by majorID_sub_big
freqy <- freqPlots(seu.obj.sub, method = 1, nrow= 1, groupBy = "majorID_sub_big", legTitle = "Cell source",refVal = "name2", showPval = F,
                   namez = unique(seu.obj$name2), 
                   colz = unique(seu.obj$colz)
                  )

freqy + ggpubr::stat_compare_means(method = "t.test",
                                   method.args = list(var.equal = F),
                                   aes(label = paste0("p = ", ..p.format..)), label.x.npc = "left", label.y.npc = 1,vjust = -1, size = 3)


ggsave(paste("./output/", outName, "/",outName, "_freqPlots.png", sep = ""), width = 10, height = 5)


#confirm variability is consistent between CIE and H
res.ftest <- lapply(levels(freqy$data$majorID_sub_big), function (x){
    data.df <- freqy$data[freqy$data$majorID_sub_big ==  x, ]
    car::leveneTest(freq ~ cellSource, data = data.df)
})


### Fig extra - DEG between Trm and Tinf
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "majorID_sub_big", idents.1 = "Tissue resident", idents.2 = "Non-resident", bioRep = "name2",padj_cutoff = 0.05, lfcCut = 0.58, 
                      minCells = 25, outDir = paste0("./output/", outName, "/"), 
                      title = "Tres_vs_Tnonres", idents.1_NAME = "Tres", idents.2_NAME = "Tnonres", 
                      returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 24
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in Tres", leftLab = "Up in Tnonres", arrowz = T) + labs(x = "log2(FC) Tres vs Tnonres") + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_TresVTnonres.png", sep = ""), width = 7, height = 7)



### Fig 2d - DEG between Trm and Tinf
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_final", idents.1 = "0", idents.2 = "1", bioRep = "name2",padj_cutoff = 0.05, lfcCut = 0.58, 
                      minCells = 25, outDir = paste0("./output/", outName, "/"), 
                      title = "CD8-TRM_vs_IL7R_T", idents.1_NAME = "Trm", idents.2_NAME = "Tinf", 
                      returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 24
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in CD8_TRM (c0)", leftLab = "Up in IL7R+_T (c1)", arrowz = T) + labs(x = "log2(FC) CD8_TRM vs IL7R+_T") + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_fig2d.png", sep = ""), width = 7, height = 7)


### Fig 2e: gsea of the DGE results
p <- plotGSEA(pwdTOgeneList = "./output/tcell/Trm_vs_Tinf_all_genes.csv", category = "C5", subcategory = NULL, 
              upCol = "blue", dwnCol = "red", size = 4.5, Trim) 

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
             label = "Up in IL7R+_T (c1)",
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
             label = "Up in CD8_TRM (c0)",
             hjust = 0.5,
             vjust = 1.5,
             size = 5)


ggsave(paste("./output/", outName, "/", outName, "_fig2e.png", sep = ""), width = 10, height = 7)



### Fig extra: gsea of the DGE results using majorID_sub_big
p <- plotGSEA(pwdTOgeneList = "./output/tcell/Tres_vs_Tnonres_all_genes.csv", category = "C7", subcategory = NULL, 
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


ggsave(paste("./output/", outName, "/", outName, "_gseaTresVTnonres.png", sep = ""), width = 10, height = 7)


########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end T cell analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

### Fig extra: dge analysis all cells
seu.obj$allCells <- "DGE analysis of all cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
linDEG(seu.obj = seu.obj, threshold = 1, thresLine = F, groupBy = "allCells", comparision = "cellSource", contrast = c("CIE", "Healthy"),
       outDir = paste0("./output/", outName,"/"), 
       outName = "all_cells", cluster = NULL, labCutoff = 15, noTitle = F,
                   colUp = "red", colDwn = "blue", subtitle = T, returnUpList = F, returnDwnList = F, forceReturn = F, useLineThreshold = F, pValCutoff = 0.01, saveGeneList = T, addLabs = ""
                  )


### Fig extra: dge analysis by major cell types
linDEG(seu.obj = seu.obj, groupBy = "majorID_sub_big", comparision = "cellSource", contrast= c("CIE","Healthy"),
       outDir = paste0("./output/", outName,"/"), outName = "majorID_sub_big", 
       pValCutoff = 0.01, saveGeneList = T, addLabs = "", labsHide = "^ENSCAFG"
                  )
