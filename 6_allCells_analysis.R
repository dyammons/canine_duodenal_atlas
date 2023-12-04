#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")

### complete analysis with n=4 cie

### Begin analysis
seu.obj <- readRDS("./output/s3/230706_duod_h3c6_NoIntrons_res0.4_dims50_dist0.5_neigh40_S3.rds")

#load in meta data
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident", metaAdd = "exclude")
seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy", "CIE"))
outName <- "allCells"
colArray <- read.csv("./refColz.csv", header = T)

table(seu.obj$exclude)
seu.obj <- subset(seu.obj,
                  subset = 
                  exclude ==  "keep")

seu.obj <- indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = "230706_duod_h3c4_NoIntrons", preSub = T, nfeatures = 2500,
                      vars.to.regress = "percent.mt"
                       )

seu.obj <- readRDS(file = "./output/s2/230706_duod_h3c4_NoIntrons_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "230706_duod_h3c4_NoIntrons", test_dims = 40, algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230706_duod_h3c4_NoIntrons", final.dims = 40, final.res = 0.3, stashID = "clusterID2", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 50, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_QC_feats.png", sep = ""), width = 9, height = 3)



#note: cluster 10 looks susspect -- filtering out (low nFeature and no clear defining feats)
seu.obj <- subset(seu.obj, invert = T,
                  subset = 
                  clusterID2 ==  "10")
table(seu.obj$clusterID2)


seu.obj <- indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = "230816_duod_h3c4_NoIntrons", preSub = T, nfeatures = 2500,
                      vars.to.regress = "percent.mt"
                       )

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "230816_duod_h3c4_NoIntrons", final.dims = 40, final.res = 1.3, stashID = "clusterID_2_1", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 50, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


#################################################
##### Analysis of n3 healthy & n4 CIE douds #####
#################################################

#load in processed data
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

seu.obj <- readRDS("./output/s3/canine_duodenum_annotated.rds")

### Fig : Create UMAP by majorID_sub
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
p <- formatUMAP(plot = pi) + NoLegend() + theme(axis.title = element_blank(),
                                                                                  panel.border = element_blank())
ggsave(paste("./output/", outName,"/",outName, "_ctl3_UMAP.png", sep = ""), width = 7, height = 7)


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
ggsave(paste("./output/", outName, "/", outName, "_ctl3_UMAP_tiny.png", sep = ""), width = 7, height = 7)



# seu.obj.all <- seu.obj

# tableData <- as.data.frame(table(seu.obj.all$celltype.l3, seu.obj.all$orig.ident))

# statData <- tableData %>% group_by(Var2) %>% mutate(samSize = sum(Freq),
#                                                    pct = Freq/samSize*100,
#                                                    cellSource = ifelse(grepl("norm",Var2),"Healthy","CIE")) %>% group_by(Var1) %>% as.data.frame()


# tstRes <- ggpubr::compare_means(pct ~ cellSource, statData, group.by = "Var1") %>% as.data.frame()



# coorData <- tableData %>% group_by(Var2) %>% mutate(samSize = sum(Freq),
#                                                    pct = Freq/samSize*100,
#                                                    cellSource = ifelse(grepl("norm",Var2),"Healthy","CIE")) %>% as.data.frame()


# allData <- statData %>% group_by(Var1) %>% summarize(avg = mean(pct),
#                                   std = sd(pct)
#                                  ) %>% mutate(cellSource = "All") %>% relocate(cellSource, .after = Var1) %>% as.data.frame()

# groupData <- read.csv("tableGroups.csv")
# headerData_all_stat <- tableData %>% left_join(groupData, by = "Var1") %>% group_by(tableGroup,Var2) %>%  summarise(Freq = sum(Freq)) %>% group_by(Var2) %>% mutate(samSize = sum(Freq),
#                                                    pct = Freq/samSize*100)

# headerData_all <- headerData_all_stat %>% group_by(tableGroup) %>% summarize(avg = mean(pct),
#                                                                              std = sd(pct)
#                                                                             ) %>% mutate(cellSource = "All") %>% relocate(cellSource, .after = tableGroup) %>% as.data.frame()

# tableData$cellSource <- ifelse(grepl("norm", tableData$Var2), "Healthy", "CIE")
# subData <- tableData %>% group_by(Var2) %>% mutate(samSize = sum(Freq),
#                                                    pct = Freq/samSize*100) %>% group_by(Var1, cellSource) %>% summarize(avg = mean(pct),
#                                                                                                             std = sd(pct)
#                                                                                                            ) %>% as.data.frame()
        
        
# headerData_sub_stat <- tableData %>% left_join(groupData, by = "Var1") %>% group_by(tableGroup,Var2) %>%  summarise(Freq = sum(Freq)) %>% group_by(Var2) %>% mutate(samSize = sum(Freq),
#                                                    pct = Freq/samSize*100,
#                                                    cellSource = ifelse(grepl("norm", Var2), "Healthy", "CIE"))
        
# headerData_sub <- headerData_sub_stat %>% group_by(tableGroup, cellSource) %>% summarize(avg = mean(pct),
#                                                                                          std = sd(pct)
#                                                                                         ) 


# compiledData <- rbind(headerData_sub,headerData_all) %>% arrange(tableGroup) %>% as.data.frame()

# statData <- rbind(headerData_sub_stat,headerData_all_stat) %>% arrange(tableGroup) %>% as.data.frame()

# tstRes <- ggpubr::compare_means(pct ~ cellSource, statData, group.by = "tableGroup") %>% as.data.frame()

# compiledData <- reshape(compiledData, idvar = "tableGroup", timevar = "cellSource", direction = "wide") %>% remove_rownames() %>%
#         column_to_rownames(var="tableGroup") %>% round(.,1) %>% 
#         transmute(All = paste(avg.All," ± ",std.All,sep=""),
#                   Healthy = paste(avg.Healthy," ± ",std.Healthy,sep=""),
#                   CIE = paste(avg.CIE," ± ",std.CIE,sep="")
#               )

# compiledData <- compiledData %>% rownames_to_column() %>% left_join(tstRes[,c("tableGroup","p.format")], by = c("rowname" = "tableGroup")) %>% remove_rownames %>% column_to_rownames(var="rowname")
        
# write.csv(compiledData, file = paste0("./output/", outName, "/", "header_compiledTableData.csv"))





# compiledData <- rbind(subData,allData) %>% arrange(Var1)
        
# compiledData <- reshape(compiledData, idvar = "Var1", timevar = "cellSource", direction = "wide") %>% filter(avg.Healthy != 0) %>% 
#         column_to_rownames(var="Var1") %>% round(.,1) %>% 
#         transmute(All = paste(avg.All," ± ",std.All,sep=""),
#                   Healthy = paste(avg.Healthy," ± ",std.Healthy,sep=""),
#                   CIE = paste(avg.CIE," ± ",std.CIE,sep="")
#               )

# compiledData <- compiledData %>% rownames_to_column() %>% left_join(tstRes[,c("Var1","p.format")], by = c("rowname" = "Var1")) %>% remove_rownames %>% column_to_rownames(var="rowname")

# write.csv(compiledData, file = paste("./output/", outName, "/", "compiledTableData.csv", sep = ""))


######################################################
### MAKE FINAL TABLE | EXTRACT MORE DATA FOR TABLE ###
######################################################

# compiledData <- rbind(headerData_sub,headerData_all) %>% arrange(majorGroup) %>% as.data.frame()

# statData <- rbind(headerData_sub_stat,headerData_all_stat) %>% arrange(majorGroup) %>% as.data.frame()

# tstRes <- ggpubr::compare_means(pct ~ cellSource, statData, group.by = "majorGroup") %>% as.data.frame()

# compiledData <- reshape(compiledData, idvar = "majorGroup", timevar = "cellSource", direction = "wide") %>% remove_rownames() %>%
#         column_to_rownames(var="majorGroup") %>% round(.,1) %>% 
#         transmute(All = paste(avg.All," ± ",std.All,sep=""),
#                   Healthy = paste(avg.Healthy," ± ",std.Healthy,sep=""),
#                   Osteosarcoma = paste(avg.Osteosarcoma," ± ",std.Osteosarcoma,sep="")
#               )

# compiledData <- compiledData %>% rownames_to_column() %>% left_join(tstRes[,c("majorGroup","p.format")], by = c("rowname" = "majorGroup")) %>% remove_rownames %>% column_to_rownames(var="rowname")
        
# write.csv(compiledData, file = paste("./output/", outName, "/", outName, "header_compiledTableData.csv", sep = ""))



#resume

#################################################
##### Analysis of n3 healthy & n4 CIE douds #####
#################################################
#ran irc using 230713_duod_h3c4_NoIntrons_res0.4_dims40_dist0.3_neigh50_S3.rds
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


saveName <- "231022_allCells_duod_h3c4_NoIntrons_2500"
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_pertyName", numOfFeats = 24, outName = saveName,
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )


ExportToCB_cus(seu.obj = seu.obj, dataset.name = outName, outDir = "./output/cb_input/", 
               markers = paste0("./output/viln/",outName,"/",saveName,"_gene_list.csv"), 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", "majorID",
                                                   "clusterID_sub", "name2", "majorID_pertyName", "cellSource", "clusterID_final"), 
               skipEXPR = F,test = F,
               feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                         "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                         "CD4", "MS4A1", "RELN","EPCAM")
                          
                          )

# Idents(seu.obj) <- "name"
# seu.obj <- RenameIdents(seu.obj, c("CIE_2" = "CIE_1", "CIE_3" = "CIE_2", 
#                                    "CIE_4" = "CIE_3")
#                        )


# seu.obj$name <- Idents(seu.obj)

# Idents(seu.obj) <- "finalHighClusID"
# seu.obj <- RenameIdents(seu.obj, c("10" = "9", "11" = "10", "12" = "11",
#                                    "13" = "12", "14" = "13",
#                                    "15" = "14","16" = "15",
#                                    "17" = "16", "18" = "17",
#                                    "19" = "18", "20" = "19"
#                                   )
#                        )


# seu.obj$finalHighClusID <- Idents(seu.obj)

# #fix a mistake
# # seu.obj$finalHighClusID <- ifelse(seu.obj$clusterID_2_1 != "19", as.character(seu.obj$finalHighClusID), "20")



# seu.obj$finalHighClusID <- factor(seu.obj$finalHighClusID, levels = as.character(sort(as.numeric(unique(seu.obj$finalHighClusID)))))


# seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./colorID_cie3.csv", groupBy = "clusterID2", metaAdd = "color")
# seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./colorID_cie3.csv", groupBy = "clusterID2", metaAdd = "labCol")
# seu.obj$cellSource <- factor(seu.obj$cellSource, levels = c("Healthy","CIE"))


#check QC params
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


#generate violin plots for each cluster
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_2_1", numOfFeats = 24, outName = "230816_duod_h3c4_NoIntrons",
                     outDir = "./output/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )


#generate violin plots for each cluster
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_pertyName", numOfFeats = 24, outName = "230816_duod_h3c4_NoIntrons",
                     outDir = "./output/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )



#plot inital cluster umap
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_final",
              cols = c("#FF89B3", "#C89504", "#00ABFD", "#B983FF", "#FA7476", "#0A9B9F", "#9B8EFF"),
              pt.size = 0.25,
              label = T,
              label.box = T
 )
p <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.8, labCol = "black") + NoLegend()+ theme(axis.title = element_blank(),
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
ggsave(paste("./output/", outName,  "/", subname,"/", outName, "_raw_UMAP_tiny.png", sep = ""), width = 7, height = 7)



### Use conocial markers to ID cells
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

# c("#FF89B3", "#C89504", "#00ABFD", "#B983FF", "#FA7476", "#0A9B9F", "#9B8EFF"
# c("T cell", "Epithelial", "Myeloid", "Plasma cell", "Cycling T cell", "Mast cell", "B cell")
colz <- c("black", "#FF89B3", "#FF89B3", 
          "#C89504", "#C89504", "#C89504",
          "#C89504", "#00ABFD", "#00ABFD",
          "#B983FF", "#FA7476", "#FA7476", 
          "#0A9B9F", "#9B8EFF", "#9B8EFF"
         )

p_feats <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, features = features, color = colz, order = T, titles = titlez, returnPlots = F, title.size = 16, pt.size = 0.00000001, min.cutoff = "q10", noLegend = T) 
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_featPlots.png", sep = ""), width = 9, height = 15)



#plot features
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
ggsave(paste("./output/", outName, "/", subname, "/", outName, "_majorDot.png", sep = ""), width = 8, height = 4)



Idents(seu.obj) <- "cellSource"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data$cellSource)))

seu.obj.ds$majorID_pertyName <- factor(seu.obj.ds$majorID_pertyName, levels = rev(c("Epithelial","T cell","Myeloid","Plasma cell","B cell","Mast cell","Cycling T cell")))
### Fig xx: stacked bar graph by colorID
p <- stackedBar(seu.obj = seu.obj.ds, downSampleBy = "cellSource", groupBy = "name2", clusters = "majorID_pertyName") +
scale_fill_manual(labels = levels(seu.obj$name2), 
                  values = levels(seu.obj$colz)) + theme(axis.title.y = element_blank(),
                                                         axis.title.x = element_text(size = 14),
                                                         axis.text = element_text(size = 12)
                                                        ) 
ggsave(paste("./output/", outName,"/", subname, "/",outName, "_stackedBar.png", sep = ""), width =7, height = 4)



### Fig 1e: stats by cie vs healthy
freqy <- freqPlots(seu.obj, method = 2, nrow= 1, groupBy = "majorID_pertyName", legTitle = "Cell source",refVal = "name2", showPval = T,
              namez = "name2", 
              colz = "colz"
              )

ggsave(paste("./output/", outName, "/", subname, "/",outName, "_freqPlots.png", sep = ""), width = 8.5, height = 3)



res.ftest <- lapply(levels(freqy$data$majorID_pertyName), function (x){
    data.df <- freqy$data[freqy$data$majorID_pertyName ==  x, ]
    leveneTest(freq ~ cellSource, data = data.df)
})
                    


### Fig 1d: umap by sample
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
fig1d <- formatUMAP(pi) + labs(colour="") + theme(legend.position = "top", legend.direction = "horizontal",legend.title=element_text(size=12)) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))
ggsave(paste("./output/", outName, "/", subname, "/", outName, "_umap_bySample.png", sep = ""), width =7, height = 7)






linDEG(seu.obj = seu.obj, groupBy = "majorID_pertyName", comparision = "cellSource", contrast= c("CIE","Healthy"),
       outDir = paste0("./output/", outName,"/", subname,"/"), outName = "majorID_pertyName", 
       pValCutoff = 0.01, saveGeneList = T, addLabs = "", labsHide = "^ENSCAFG"
                  )


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


### load in kallisto output
outName <- "kallisto"
seu.kal <- MuDataSeurat::ReadH5AD("/scratch/alpine/dyammons@colostate.edu/proj03_k9_duod/03_output/2023-11-06_kb-count_output/duod_norm_1/counts_unfiltered/adata.h5ad")

#visualize QC metrics as a violin plot
seu.kal[["percent.mt"]] <- PercentageFeatureSet(seu.kal, pattern = "^MT-")
p <- VlnPlot(seu.kal, features = c("nFeature_RNA", "nCount_RNA"))
ggsave(paste0("./output/", outName,"/", outName, "_QC.png") )

seu.obj <- readRDS("./output/s3/canine_duodenum_annotated.rds")
seu.obj <- subset(seu.obj, subset = orig.ident == "run_count_duod_norm_1_noIntrons")
goodCodez <- substr(colnames(seu.obj), 1, nchar(colnames(seu.obj))-4)

p <- VlnPlot(seu.obj, features = c("nFeature_RNA", "nCount_RNA"))
ggsave(paste0("./output/", outName,"/", outName, "_QC_orig.png") )


unname(seu.kal$barcode) %in% goodCodez %>% table()
#  FALSE   TRUE 
# 662842   4745 

seu.kal <- subset(seu.kal, subset = barcode %in% goodCodez)
p <- VlnPlot(seu.kal, features = c("nFeature_RNA", "nCount_RNA"))
ggsave(paste0("./output/", outName,"/", outName, "_QC.png") )

# ### Do not use this code block
# #in terminal -- cut -f1-3 transcripts_to_genes.txt | sed 's/\t/_/g' > convert.txt
# tr2g.df <- read.table("/projects/dyammons@colostate.edu/references/canine/convert.txt", header = F)
# tr2g.df <- tr2g.df %>% separate(V1, c('trans', 'ens','gene_symbol'), sep = "_") %>% mutate(gene_symbol = ifelse(gene_symbol == "", ens, gene_symbol)) 
# tr2g.df <- tr2g.df[ ,-1] %>% distinct()


# genez <- as.data.frame(rownames(seu.kal))
# colnames(genez) <- "ens"
# newGenez <- genez %>% left_join(tr2g.df, by = "ens")

# gene symbols are stashed in meta.features slot, so extract them
newGenez <- seu.kal@assays$RNA@meta.features %>% rownames_to_column("gene") %>% mutate(newGenez = ifelse(duplicated(seu.kal@assays$RNA@meta.features$gene_name), as.character(gene), as.character(gene_name))) %>% pull(newGenez)

cnt.mat <- as.matrix(seu.kal@assays$RNA@counts)
colnames(cnt.mat) <- goodCodez
rownames(cnt.mat) <- newGenez

meta.df.orig <- seu.obj@meta.data
meta.df.orig$barcode <- goodCodez

#test that the order matches up
table(rownames(meta.df.orig) == paste0(meta.df.orig$barcode,"-1_7"))
# TRUE 
# 4745 

meta.df <- seu.kal@meta.data %>% left_join(meta.df.orig, by = "barcode")
rownames(meta.df) <- meta.df$barcode
# Generate new Seurat object.
seu.obj <- CreateSeuratObject(
     cnt.mat,
     project = "SeuratProject",
     assay = "RNA",
     min.cells = 0,
     min.features = 0,
     names.field = 1,
     names.delim = "_",
     meta.data = meta.df
)


# reduce the data
seu.obj <- NormalizeData(seu.obj,
                         normalization.method = "LogNormalize",
                         Scale.factor = 10000) ###change method to scran???

seu.obj <- FindVariableFeatures(seu.obj,
                                selection.method = "vst", 
                                nfeatures = 2500) #can change number of feats used

all.genes <- rownames(seu.obj)
seu.obj <- ScaleData(seu.obj, features = all.genes)
seu.obj <- RunPCA(seu.obj, features = VariableFeatures(object = seu.obj))
seu.obj <- FindNeighbors(seu.obj,
                         dims = 1:10
                        ) #can change dims
seu.obj <- FindClusters(seu.obj,
                        resolution = 0.1
                       ) #can change resolution
seu.obj <- RunUMAP(seu.obj, 
                   dims = 1:15
                  ) #can change dims


featPlots = c("PTPRC", "CD3E", "CD8A", "GZMA", 
              "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
              "CD4", "MS4A1", "TRDC","FOXP3")

### Fig : Create UMAP by majorID_sub
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "celltype.l3",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = TRUE
)

ggsave(paste0("./output/", outName,"/",outName, "_rawUMAP.png"), width = 7, height = 7)




seu.obj <- readRDS("./output/s3/canine_duodenum_annotated.rds")
seu.obj <- subset(seu.obj, subset = cellSource == "Healthy")
vilnPlots(seu.obj = seu.obj, groupBy = "celltype.l3", numOfFeats = 24, outName = "231108_duod_h3_cr",
                     outDir = paste0("./output/viln/allCells/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )


loadKal <- function(din = NULL, dout = NULL, seu.obj = NULL){

    fpath <-  paste0("./", din,"/") 

    files <- list.files(path = fpath, pattern=NULL, all.files=FALSE,
                        full.names=F)
    
    for (infile in files) {
        pwd <- paste0("./", din,"/", infile, "/adata.h5ad") 

        seu.kal <- MuDataSeurat::ReadH5AD(pwd)

        # gene symbols are stashed in meta.features slot, so extract them
        dups <- as.character(seu.kal@assays$RNA@meta.features[duplicated(seu.kal@assays$RNA@meta.features$gene_name), ])
        if(length(dups) > 0){
            message(paste0("Duplicate gene symbols found in the seu.kal@assays$RNA@meta.features slot. A total of ", length(dups), " duplicated gene symbols were identified and the gene ID will be retained.\n"))
            message(paste0("The impacted features are: ", paste(dups, collapse = ", ")))

            newGenez <- seu.kal@assays$RNA@meta.features %>% rownames_to_column("gene") %>% mutate(newGenez = ifelse(duplicated(seu.kal@assays$RNA@meta.features$gene_name), as.character(gene), as.character(gene_name))) %>% pull(newGenez)

        }

        message(paste0("Subseting the annotated dataset for transfer.", 
                       "Searching for ",infile," in ", paste(as.character(unique(seu.obj$orig.ident)) , collapse = ", "))
               )
        seu.obj.sub <- subset(seu.obj, cells = colnames(seu.obj)[grepl(infile, seu.obj$orig.ident)])
        goodCodez <- substr(colnames(seu.obj.sub), 1, 16) #barcodes must be nchar == 16

        message("Subseting the kallisto pseudo-aligned object to remove empty droplets and other low quality data points.")
        seu.kal <- subset(seu.kal, subset = barcode %in% goodCodez)

        message(paste0("Of the ", nrow(seu.obj.sub@meta.data), " in the annotated dataset, ", nrow(seu.kal@meta.data), " barcodes were also found in the kallisto pseudo-aligned sample."))
        
        cnt.mat <- as.matrix(seu.kal@assays$RNA@counts)
        colnames(cnt.mat) <- goodCodez
        rownames(cnt.mat) <- newGenez

        meta.df.orig <- seu.obj.sub@meta.data
        meta.df.orig$barcode <- goodCodez

        #test that the order matches up
        table(substr(rownames(meta.df.orig), 1, 16) == substr(colnames(seu.obj.sub), 1, 16))
        # TRUE 
        # 4745 

        meta.df <- seu.kal@meta.data %>% left_join(meta.df.orig, by = "barcode")
        rownames(meta.df) <- meta.df$barcode

        #generate new Seurat object
        seu.obj.kal <- CreateSeuratObject(cnt.mat,
                                      project = infile,
                                      assay = "RNA",
                                      min.cells = 0,
                                      min.features = 0,
                                      names.field = 1,
                                      names.delim = "_",
                                      meta.data = meta.df
                                     )


        #reduce the data
        seu.obj.kal <- NormalizeData(seu.obj.kal,
                                 normalization.method = "LogNormalize",
                                 Scale.factor = 10000) ###change method to scran???

        seu.obj.kal <- FindVariableFeatures(seu.obj.kal,
                                        selection.method = "vst", 
                                        nfeatures = 2500) #can change number of feats used

        all.genes <- rownames(seu.obj.kal)
        seu.obj.kal <- ScaleData(seu.obj.kal, features = all.genes)
        seu.obj.kal <- RunPCA(seu.obj.kal, features = VariableFeatures(object = seu.obj.kal))
        seu.obj.kal <- FindNeighbors(seu.obj.kal,
                                 dims = 1:10
                                ) #can change dims
        seu.obj.kal <- FindClusters(seu.obj.kal,
                                resolution = 0.1
                               ) #can change resolution
        seu.obj.kal <- RunUMAP(seu.obj.kal, 
                           dims = 1:15
                          ) #can change dims


        featPlots = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                      "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                      "CD4", "MS4A1", "TRDC","FOXP3")

        ### Fig: Create UMAP by celltype.l3
        pi <- DimPlot(seu.obj.kal, 
                      reduction = "umap", 
                      group.by = "celltype.l3",
                      pt.size = 0.25,
                      label = T,
                      label.box = T,
                      shuffle = TRUE
        )

        ggsave(paste0("./",dout,"/",infile, "_rawUMAP.png"), width = 7, height = 7)

        saveRDS(seu.obj.kal, file = paste0("./",dout,"/", infile,"_kallisto_S1.rds"))
    }

}


### diverge
outName <- "231110_duod_cie3_cr"
seu.obj <- sctIntegrate(din = "./output/s1/", dout = "./output/s2/", outName = "231110_duod_cie3_cr", vars.to.regress = "percent.mt", nfeatures = 2500, pattern = ".*cie.*_S1.rds")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "231110_duod_cie3_cr", final.dims = 40, final.res = 0.8, stashID = "clusterID_kal", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 50, assay = "integrated", saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

seu.obj <- sctIntegrate(din = "./output/s1/", dout = "./output/s2/", outName = "231110_duod_h3_cr", vars.to.regress = "percent.mt", nfeatures = 2500, pattern = ".*norm.*_S1.rds")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "231110_duod_h3_cr", final.dims = 40, final.res = 0.8, stashID = "clusterID_kal", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 50, assay = "integrated", saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

seu.obj <- readRDS("./output/s3/231110_duod_h3_cr_res0.8_dims40_dist0.3_neigh50_S3.rds")
seu.obj.all <- readRDS("./output/s3/canine_duodenum_annotated.rds")

seu.obj.sub <- indReClus(seu.obj = seu.obj.all, outDir = "./output/s2/", subName = "231110_all_h3c4_normVar_NoIntrons_2500", 
                     preSub = T, nfeatures = 2500, vars.to.regress = "percent.mt", saveRDS = F, varFeatsToUse = seu.obj@assays$integrated@var.features
                    )

seu.obj.sub <- dataVisUMAP(seu.obj = seu.obj.sub, outDir = "./output/s3/", outName = "231110_all_h3c4_normVar_NoIntrons_2500", final.dims = 40, final.res = 0.8, stashID = "clusterID_kal", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 50, assay = "integrated", saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


seu.obj.all <- readRDS("./output/s3/canine_duodenum_annotated.rds")
seu.obj.all <- subset(seu.obj.all, subset = cellSource == "Healthy")

seu.obj.sub <- indReClus(seu.obj = seu.obj.all, outDir = "./output/s2/", subName = "231110_all_h3_NoIntrons_2500", 
                     preSub = T, nfeatures = 2500, vars.to.regress = "percent.mt", saveRDS = F,
                    )





seu.obj <- readRDS("./output/s3/canine_duodenum_annotated.rds")
loadKal(din = "./kal_Input", dout = "./output/kal_s1/",
                    seu.obj = seu.obj)



seu.obj <- sctIntegrate(din = "./output/kal_s1/", dout = "./output/kal_s2/", outName = "231108_duod_h3_cie4_kal", vars.to.regress = "percent.mt", nfeatures = 2500)


seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "231108_duod_h3_cie4_kal", final.dims = 40, final.res = 0.8, stashID = "clusterID_kal", algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 50, assay = "integrated", saveRDS = T, return_obj = T, returnFeats = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


#generate violin plots for each cluster
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_kal", numOfFeats = 24, outName = "231108_duod_h3_kal",
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )

vilnPlots(seu.obj = seu.obj, groupBy = "celltype.l3", numOfFeats = 24, outName = "231108_duod_h3_kal_ctl3",
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )


### Fig : Create UMAP by majorID_sub
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "celltype.l3",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = TRUE
)

ggsave(paste0("./output/",outName,"/", outName, "_rawUMAP.png"), width = 10, height = 7)



################################################
#                                              #
#                                              #
#         BEGIN DECONVOLUTION ANALYSIS         #
#                                              #
#                                              #
################################################


### Load in processed object
seu.obj <- readRDS("./output/s3/231108_duod_h3_kal_res0.8_dims40_dist0.3_neigh50.rds")
seu.obj <- readRDS("./output/s3/231108_duod_h3_cie4_kal_res0.8_dims40_dist0.3_neigh50_S3.rds")
outName <- "kallisto"


### Create pseudo bulk matirx to test out the decon results
seu.obj$allCells <- "allCells"
createPB(seu.obj = seu.obj, groupBy = "allCells", comp = "cellSource", biologicalRep = "name2", cnts = F, dwnSam =F,
                     clusters = NULL, outDir = "./output/kallisto/",
                     grepTerm = "H", grepLabel = c("healthy","healthy"))
                    
write.table(file = "./output/kallisto/src/data/groundTruth.txt", read.csv("./output/kallisto/allCells_pb_matrix.csv"), row.names = F, quote=T, sep='\t')


### Set level for deconvolution
#extract data to modify levels to create a higher level called "ctl2"
# write.csv(unique(seu.obj$celltype.l3), file = "./cellType_levels.csv", row.names = F)

#read in metadata and stash idents in ctl2
ct.df <- read.csv("./cellType_levels.csv")
namez <- gsub(" ", "_", ct.df$ctl2)
names(namez) <- ct.df$ctl3

Idents(seu.obj) <- "celltype.l3"
seu.obj <- RenameIdents(seu.obj, namez
                       )

seu.obj$ctl2 <- Idents(seu.obj)


namez <- NULL
namez <- gsub(" ", "_", ct.df$ctl2_2)
names(namez) <- ct.df$ctl3

Idents(seu.obj) <- "celltype.l3"
seu.obj <- RenameIdents(seu.obj, namez
                       )

seu.obj$ctl2_deCon <- Idents(seu.obj)


#generate defining features for downstream use in filtering of deconvolution matrix
#for ctl2
vilnPlots(seu.obj = seu.obj, groupBy = "ctl2", numOfFeats = 24, outName = "231116_duod_h3_kal",
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )

#for ctl2_2
#remove IFN-sig cell type to reduce skewing of data toward these populations
seu.obj.sub <- subset(seu.obj, invert = T, subset = ctl2_deCon == "remove")
seu.obj.sub$ctl2_deCon <- droplevels(seu.obj.sub$ctl2_deCon)
vilnPlots(seu.obj = seu.obj.sub, groupBy = "ctl2_deCon", numOfFeats = 24, outName = "231121_duod_h3_kal",
                     outDir = paste0("./output/viln/",outName,"/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )


### Get the true cell type percentages as determined using scRNA seq
#for ctl2
groupByList <- seu.obj$name2
clusterList <- seu.obj$ctl2

cluster_freq.table <- as.data.frame(table(groupByList, clusterList)) %>% melt()
cluster_freq.table <- cluster_freq.table[,-3]
colnames(cluster_freq.table) <- c("Sample", "ClusterID", "Count")

cluster_freq.table <- cluster_freq.table %>% dplyr::group_by(Sample) %>% mutate(pct = round(prop.table(Count),2))

#for ctl2_2
groupByList <- seu.obj.sub$name2
clusterList <- seu.obj.sub$ctl2_deCon

cluster_freq.table <- as.data.frame(table(groupByList, clusterList)) %>% melt()
cluster_freq.table <- cluster_freq.table[,-3]
colnames(cluster_freq.table) <- c("Sample", "ClusterID", "Count")

cluster_freq.table <- cluster_freq.table %>% dplyr::group_by(Sample) %>% mutate(pct = round(prop.table(Count),2))


### Export matrix to generate a reference for deconvolution using cibersortx
#for ctl2
seu.obj.sub <- subset(seu.obj, downsample = min(table(seu.obj$ctl2)))
mat <- as.matrix(GetAssayData(object = seu.obj.sub, slot = "counts"))
colnames(mat) <- gsub(" ", "_", unname(seu.obj.sub$ctl2))
mat <- cbind(gene = rownames(mat), as.data.frame(mat))

mat <- mat[rowSums(mat > 0) >= 10, ]

write.table(mat,'./output/kallisto/src/data/counts.txt', sep = '\t', row.names = F, col.names = T, quote = F)

#for ctl2_2
seu.obj.sub <- subset(seu.obj.sub, downsample = min(table(seu.obj.sub$ctl2_deCon)))
mat <- as.matrix(GetAssayData(object = seu.obj.sub, slot = "counts"))
mat <- cbind(gene = rownames(mat), as.data.frame(mat))

# #for cibdersort ref
# mat <- as.matrix(GetAssayData(object = seu.obj.sub, slot = "counts"))
# colnames(mat) <- gsub(" ", "_", unname(seu.obj.sub$ctl2_deCon))
# mat <- cbind(gene = rownames(mat), as.data.frame(mat))

# write.table(mat,'./output/kallisto/src/data/counts_ctl2_2.txt', sep = '\t', row.names = F, col.names = T, quote = F)

meta.df <- seu.obj.sub@meta.data[ ,c("name2", "ctl2_deCon")]
colnames(meta.df) <- c("sample.id", "cellType")

### Run cibersort default using Singularity
# singularity exec -B /pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/kallisto/src/data:/src/data -B /pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/kallisto/src/outdir:/src/outdir /projects/$USER/software/cibersort/fractions_latest.sif /src/CIBERSORTxFractions \
# --username dyammons@colostate.edu \
# --token e89f071b7e3ff9302a668e0f04082625 \
# --single_cell TRUE \
# --refsample /src/data/counts.txt \
# --mixture /src/data/groundTruth.txt \
# --fraction 0 \
# --rmbatchSmode TRUE 

# singularity exec -B /pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/kallisto/src/data:/src/data -B /pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/kallisto/src/outdir:/src/outdir /projects/$USER/software/cibersort/fractions_latest.sif /src/CIBERSORTxFractions \
# --username dyammons@colostate.edu \
# --token e89f071b7e3ff9302a668e0f04082625 \
# --single_cell TRUE \
# --refsample /src/data/counts.txt \
# --mixture /src/data/groundTruth.txt \
# --fraction 0.8 \
# --rmbatchSmode TRUE 

# singularity exec -B /pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/kallisto/src/data:/src/data -B /pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/kallisto/src/outdir:/src/outdir /projects/$USER/software/cibersort/fractions_latest.sif /src/CIBERSORTxFractions \
# --username dyammons@colostate.edu \
# --token e89f071b7e3ff9302a668e0f04082625 \
# --single_cell TRUE \
# --refsample /src/data/counts.txt \
# --mixture /src/data/groundTruth.txt \
# --fraction 0.8 \
# --rmbatchSmode FALSE 


# singularity exec -B /pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/kallisto/src/data:/src/data -B /pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/kallisto/src/outdir:/src/outdir /projects/$USER/software/cibersort/fractions_latest.sif /src/CIBERSORTxFractions \
# --username dyammons@colostate.edu \
# --token e89f071b7e3ff9302a668e0f04082625 \
# --single_cell TRUE \
# --refsample /src/data/counts.txt \
# --mixture /src/data/groundTruth.txt \
# --G.min 100 \
# --G.max 300 \
# --fraction 0.8 \
# --rmbatchSmode FALSE 

# singularity exec -B /pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/kallisto/src/data:/src/data -B /pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/kallisto/src/outdir:/src/outdir /projects/$USER/software/cibersort/fractions_latest.sif /src/CIBERSORTxFractions \
# --username dyammons@colostate.edu \
# --token e89f071b7e3ff9302a668e0f04082625 \
# --single_cell TRUE \
# --refsample /src/data/counts.txt \
# --mixture /src/data/groundTruth.txt \
# --G.min 50 \
# --G.max 100 \
# --fraction 0.8 \
# --rmbatchSmode FALSE


# #what worked well for decon of artifical bulk data does not work well for true bulk data :|

# singularity exec -B /pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/kallisto/src/data:/src/data -B /pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/kallisto/src/outdir:/src/outdir /projects/$USER/software/cibersort/fractions_latest.sif /src/CIBERSORTxFractions \
# --username dyammons@colostate.edu \
# --token e89f071b7e3ff9302a668e0f04082625 \
# --single_cell TRUE \
# --refsample /src/data/counts.txt \
# --mixture /src/data/groundTruth.txt \
# --G.min 50 \
# --G.max 150 \
# --fraction 0.25 \
# --rmbatchSmode TRUE \
# --remake TRUE \
# --replicates 50



#see script

### Evaluate cibersortX performance

#load in cibersort output
frac <- c("frac_0","frac_0.10","frac_0.20","frac_0.30","frac_0.40","frac_0.50","frac_0.60","frac_0.70","frac_0.80","frac_0.90","frac_1")
frac <- "ctl2_2"

filez <- paste0("./output/kallisto/", frac, "/src/outdir/CIBERSORTx_Adjusted.txt")
# filez <- paste0("./output/kallisto/", frac, "/src/outdir/CIBERSORTx_Results.txt")


cbr_res_df_list <- lapply(filez, read.table, header = T)

cbr_res_df <- do.call(rbind, cbr_res_df_list)

cbr_res_df$frac <- rep(frac, each = nrow(cbr_res_df_list[[1]]))


res <- lapply(1:length(frac), function(x){
    
    #extract one fraction value
    cbr_res_df <- cbr_res_df_list[[x]]
    saveName <- frac[x]
    
    #get the data together to plot
    df_toPlot <- cbr_res_df[ ,1:(ncol(cbr_res_df)-3)] %>% melt() %>% left_join(cluster_freq.table, by = c("Mixture" = "Sample", "variable" = "ClusterID"))

    df_toPlot <- df_toPlot[!is.na(df_toPlot$pct), ]

    rmse <- sqrt(mean((df_toPlot$pct - df_toPlot$value)^2)) * 100

    # df_toPlot %>% group_by(variable) %>% summarize(RMSE = sqrt(mean((pct - value)^2))*100)
    # df_toPlot %>% group_by(variable) %>% summarize(RMSE = sqrt(mean((abs((pct-value)/(pct+1)) - pct)^2)))


#     df_fit <- df_toPlot %>% group_by(variable) %>%
#       do(fit = lm(value ~ pct, data = .))


#     df_toPlot %>% nest_by(variable) %>%
#       mutate(mod = list(lm(value ~ pct, data = df_toPlot))) %>%
#       dplyr::reframe(broom::tidy(mod))


#     df_fit_Coef = broom::tidy(df_fit, fit)
#     df_fit_Coef

    #plot the data
    p <- ggplot(df_toPlot, aes(x=pct, y=value)) + 
                    stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linetype = "dashed", color = "grey50") +
                    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
                    geom_point() +
                    scale_x_continuous(limits = c(0, 1)) +
                    scale_y_continuous(limits = c(0, 1)) +
                    labs(x = "Truth (pct)", y = "Predicted (pct)") +
                    guides(color = guide_legend(title = "Cell type", size = 3, override.aes=list(fill=NA))) +
                    #geom_text(x = 75, y = 5, label = as.character(as.expression(eq)), parse = TRUE) +
                    theme(panel.background = element_rect(fill = "transparent",colour = NA),
                          plot.background = element_rect(fill = "transparent",colour = NA),
                          legend.background = element_rect(fill = "transparent",colour = NA),
                          legend.key = element_rect(fill = "transparent",colour = NA),
                          panel.grid.major = element_line(color = "gray"), 
                          panel.grid.minor = element_line(color = "gray"),#axis.text = element_blank(), 
                          axis.ticks = element_blank(),
                          axis.title = element_text(size= 14),
                          #plot.title = element_text(hjust = 0.02),
                          #plot.title.position = "plot",
                          plot.title = element_blank(),
                          title = element_text(size= 16),
                          axis.line = element_blank(),
                          panel.border = element_rect(color = "black",
                                                      fill = NA,
                                                      size = 2)
                          ) + stat_cor()

    pi <- p + facet_wrap("variable", scales = "free_y", nrow = 8) + stat_cor() + theme(plot.background = element_rect(fill = "white"))
    ggsave(paste0("./output/",outName,"/", outName, "_", saveName, "corrPlots.png"), width = 7, height = 10)
    
    return(rmse)

})








### methods do not account for transcriptome size... so use deviatoin from average to scale the data
avg <- mean(seu.obj$nFeature_RNA)
meta.df <- seu.obj@meta.data

scalez <- meta.df %>% group_by(ctl2) %>% summarize(scalez = mean(nFeature_RNA)/avg) %>% as.data.frame()


deconRES_scaled <- df_toPlot %>% left_join(scalez, by = c("variable" = "ctl2")) %>% mutate(scaledValue = value/scalez)
fixIT <- deconRES_scaled %>% group_by(Mixture) %>% summarize(chck = sum(scaledValue))
deconRES_scaled_correct <- deconRES_scaled %>% left_join(fixIT, by = "Mixture") %>% mutate(pctCorrect = (1/chck)*scaledValue)
deconRES_scaled_correct %>% group_by(Mixture) %>% summarize(chck = sum(pctCorrect))



#plot the data
p <- ggplot(deconRES_scaled_correct, aes(x=pct, y=pctCorrect)) + 
                stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linetype = "dashed", color = "grey50") +
                geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
                geom_point() +
                scale_x_continuous(limits = c(0, 1)) +
                scale_y_continuous(limits = c(0, 1)) +
                labs(x = "Truth (pct)", y = "Predicted (pct)") +
                guides(color = guide_legend(title = "Cell type", size = 3, override.aes=list(fill=NA))) +
                #geom_text(x = 75, y = 5, label = as.character(as.expression(eq)), parse = TRUE) +
                theme(panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA),
                      legend.background = element_rect(fill = "transparent",colour = NA),
                      legend.key = element_rect(fill = "transparent",colour = NA),
                      panel.grid.major = element_line(color = "gray"), 
                      panel.grid.minor = element_line(color = "gray"),#axis.text = element_blank(), 
                      axis.ticks = element_blank(),
                      axis.title = element_text(size= 14),
                      #plot.title = element_text(hjust = 0.02),
                      #plot.title.position = "plot",
                      plot.title = element_blank(),
                      title = element_text(size= 16),
                      axis.line = element_blank(),
                      panel.border = element_rect(color = "black",
                                                  fill = NA,
                                                  size = 2)
                      ) + stat_cor()
                
pi <- p + facet_wrap("variable", scales = "free_y", nrow = 8) + stat_cor() + theme(plot.background = element_rect(fill = "white"))

saveName <- "ciberBase"
ggsave(paste0("./output/",outName,"/", outName, "_", saveName, "corrPlots_scaled.png"), width = 7, height = 10)

























### Create a reference for deconvolution
createCIBERsort(seu.obj = seu.obj, groupBy = "ctl2", downSample = T, outDir = paste0("./output/", outName,"/"), outName = "kallisto_ctl2", normMethod = "tpm"
                    )

#filter the reference using viln output -- ref 1
classy.df <- read.csv("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/viln/kallisto/231116_duod_h3_kal_gene_list.csv", row.names = 1, header = T)
# classy.df <- classy.df %>% filter(duplicated(gene) == FALSE)
classy.df <- classy.df[classy.df$avg_log2FC > 0.58 & classy.df$pct.1 > 0.50, ]

ds_val <- classy.df %>% group_by(cluster) %>% summarize(cntz = n())

#load a filter .csv files for ciber sort
ciberDS.df <- read.csv("./output/kallisto/kallisto_ctl2_ciberSort_matrix.csv", row.names = 1, header = T)
ciberDS.df.filtered <- ciberDS.df[row.names(ciberDS.df) %in% classy.df$gene,]
ciberDS.df.filtered$gene <- rownames(ciberDS.df.filtered)
ciberDS.df.filtered <- ciberDS.df.filtered %>% relocate(gene)
write.table(ciberDS.df.filtered, file="./output/kallisto/src/data/kallisto_ctl2_ciberSort_matrix_ref1.txt", row.names = F, col.names = T, quote=T, sep='\t')



#filter the reference using viln output -- ref 2
classy.df <- read.csv("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/viln/kallisto/231116_duod_h3_kal_gene_list.csv", row.names = 1, header = T)
# classy.df <- classy.df %>% filter(duplicated(gene) == FALSE)
classy.df <- classy.df[classy.df$avg_log2FC > 0.58 & classy.df$pct.1 > 0.50, ]

ds_val <- classy.df %>% group_by(cluster) %>% summarize(cntz = n()) %>% pull(cntz) %>% min()

classy.df <- classy.df %>% group_by(cluster) %>% top_n(., -ds_val, wt = p_val_adj)

#load a filter .csv files for ciber sort
ciberDS.df <- read.csv("./output/kallisto/kallisto_ctl2_ciberSort_matrix.csv", row.names = 1, header = T)
ciberDS.df.filtered <- ciberDS.df[row.names(ciberDS.df) %in% classy.df$gene,]
ciberDS.df.filtered$gene <- rownames(ciberDS.df.filtered)
ciberDS.df.filtered <- ciberDS.df.filtered %>% relocate(gene)
write.table(ciberDS.df.filtered, file="./output/kallisto/src/data/kallisto_ctl2_ciberSort_matrix_ref2.txt", row.names = F, col.names = T, quote=T, sep='\t')






singularity exec -B /pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/kallisto/src/data:/src/data -B /pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/kallisto/src/outdir:/src/outdir /projects/$USER/software/cibersort/fractions_latest.sif /src/CIBERSORTxFractions \
--username dyammons@colostate.edu \
--token e89f071b7e3ff9302a668e0f04082625 \
--sigmatrix /src/data/kallisto_ctl2_ciberSort_matrix.txt \
--mixture /src/data/groundTruth.txt

































### Load in the bulk data
df <- read.csv("./bulk_decoder.csv")
df_samNames <- df %>% mutate(namez = ifelse(sc_name == "na",bulk_name, sc_name))

#https://nbisweden.github.io/workshop-RNAseq/2011/lab_kallisto.html
# #in terminal -- cut -f1-3 transcripts_to_genes.txt | sed 's/\t/_/g' > convert.txt

#### NOTE THIS FILE IS MISSING 613 transcripts
tr2g.df <- read.table("/projects/dyammons@colostate.edu/references/canine/convert.txt", header = F)
tr2g.df <- tr2g.df %>% separate(V1, c('feature', 'ens','gene_symbol'), sep = "_") %>% mutate(gene_symbol = ifelse(gene_symbol == "", ens, gene_symbol)) 

#load in & save the cnt matricies
library(tximport)

files <- paste("/scratch/alpine/dyammons@colostate.edu/proj03_k9_duod/03_output/2023-11-03_output/02_kallisto" ,  
               list.files(path = "/scratch/alpine/dyammons@colostate.edu/proj03_k9_duod/03_output/2023-11-03_output/02_kallisto", pattern = "abundance.tsv", recursive = TRUE),
               sep = "/")

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tr2g.df, ignoreAfterBar = TRUE)

files <- list.files(path = "/scratch/alpine/dyammons@colostate.edu/proj03_k9_duod/03_output/2023-11-03_output/02_kallisto", pattern=NULL, all.files=FALSE,
                        full.names=F)

cts <- as.data.frame(txi.kallisto.tsv$counts)
colnames(cts) <- df_samNames$namez
cts$ens <- rownames(cts)

# write.table(cts, file = "./output/kallisto/mixture.tsv", row.names=TRUE, sep="\t")

#remove duplicate genes and rename samples to match the scRNA nomenclature
tr2g.df <- tr2g.df %>% group_by(ens) %>% slice_sample(n = 1)

cts <- cts %>% left_join(tr2g.df, by = "ens")
duplicated(cts$gene_symbol) %>% table()
# FALSE  TRUE 
# 20240 17 

newGenez <- cts %>% mutate(newGenez = ifelse(duplicated(cts$gene_symbol), as.character(ens), as.character(gene_symbol))) %>% pull(newGenez)
length(newGenez) == nrow(cts)
#[1] TRUE
rownames(cts) <- newGenez

mix <- cbind(newGenez, cts[ ,1:11])

write.table(mix, file = "./output/kallisto/src/data/mixture.txt", row.names=FALSE, sep="\t")










### Run SQUID
library(SQUID) #need MAST and DWLS installed
library(DWLS)
devtools::source_url("https://github.com/favilaco/deconv_matching_bulk_scnRNA/blob/master/helper_functions.R?raw=TRUE")
createCIBERsort(seu.obj = seu.obj.sub, groupBy = "ctl2_deCon", downSample = T, outDir = paste0("./output/", outName,"/"), outName = "kallisto_ctl2_2", normMethod = "tpm"
                    )
ciberDS.df <- read.csv("./output/kallisto/kallisto_ctl2_2_ciberSort_matrix.csv", row.names = 1, header = T)

mat$gene <- NULL
mix$newGenez <- NULL
res <- SQUID(B = as.matrix(mix), scC = as.matrix(mat) , scMeta = meta.df, pB = as.matrix(ciberDS.df), P = NULL, LeaveOneOut = FALSE)
































### Load in ground truth data
mix_grdT <- read.table("./output/kallisto/groundTruth.tsv", header = T)
rownames(mix_grdT) <- mix_grdT$X
mix_grdT <- mix_grdT[ ,-1]


### Run deconvolution using EPIC
library(EPIC)

sig <- ciberDS.df.filtered
rownames(sig) <- sig$gene
sig <- sig[ ,-1]

# refProfiles.var = lapply(levels(seu.obj$ctl2), function(x){
#     mat <- seu.obj@assays$RNA@counts[, WhichCells(seu.obj, ident  = x)]
#     outz <- as.data.frame(matrixStats::rowSds(Matrix::as.matrix(mat)))
#     colnames(outz) <- x
    
#     return(outz)
# })

# refProfiles.var = do.call(cbind.data.frame, refProfiles.var)

sig.ref <- list(as.matrix(sig), rownames(sig))#, as.matrix(refProfiles.var),
names(sig.ref) <- c("refProfiles","sigGenes")#"refProfiles.var",
out <- EPIC(bulk = mix_grdT, reference = sig.ref)
deconRES <- out$cellFractions
deconRES <- t(deconRES) %>% melt()
deconRES$value <- as.numeric(deconRES$value)




#### the numbers aren't adding up -- there is something weird here
kal.df <- read.csv("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/viln/kallisto/231108_duod_h3_kal_ctl3_gene_list.csv")
cr.df <- read.csv("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/viln/allCells/231108_duod_h3_cr_gene_list.csv")

#calculate the raw nmber of overlaping DEGs
res <- lapply(unique(dog.df$cluster), function(x){
    
    dog.list <- dog.df[dog.df$cluster == x, ] %>% pull(gene)
    
    res_pre <- lapply(unique(human.df$cluster), function(y){
        human.list <- human.df[human.df$cluster == y, ] %>% pull(gene)
#         degList <- length(unique(c(dog.list,human.list)))
        intERsct <- length(dog.list[dog.list %in% human.list])
        names(intERsct) <- x
        
        return(intERsct)
    })
    
    res_pre <- do.call(rbind, res_pre)
    rownames(res_pre) <- unique(human.df$cluster)
    return(res_pre)
    
})

res <- lapply(unique(kal.df$cluster), function(x){
    
    kal.list <- kal.df[kal.df$cluster == x, ] %>% pull(gene)
    cr.list <- cr.df[cr.df$cluster == x, ] %>% pull(gene)
    
    #kal <- kal.list[!kal.list %in% cr.list]
    #cr <- cr.list[!cr.list %in% kal.list]
    kal <- length(kal.list)
    cr <- length(cr.list)
    intERsct <- intersect(kal.list,cr.list)

        return(list("kal" = kal,
                    "cr" = cr,
                    "intERsct" = intERsct))
})
    
    res_pre <- do.call(rbind, res_pre)
    rownames(res_pre) <- unique(human.df$cluster)
    return(res_pre)
    
})

for (x in unique(kal.df$cluster)) {
    print(x)
    
    kal.df.sub <- kal.df[kal.df$cluster == x, ]
    cr.df.sub <- cr.df[cr.df$cluster == x, ]
    
    df.list <- list(kal.df.sub,cr.df.sub)

    feats.list <- lapply(df.list, function(x){feats <- x %>% filter(avg_log2FC > 0) %>% select(X)})
#     feats.list <- lapply(df.list, function(x){feats <- x %>% filter(avg_log2FC > 0, p_val_adj < 0.01) %>% select(X)})


    library(UpSetR)
    upSet.df <- as.data.frame(c(feats.list[1][[1]]$X,feats.list[2][[1]]$X))
    colnames(upSet.df) <- "gene"

    upSet.df$`Kallisto` <- as.integer(ifelse(upSet.df$gene %in% feats.list[1][[1]]$X, 1, 0))
    upSet.df$`cellranger` <- as.integer(ifelse(upSet.df$gene %in% feats.list[2][[1]]$X, 1, 0))

    # Plot sample distance heatmap with ComplexHeatmap
    png(file = paste0("./output/", outName,"/", x, "_upSet.png"), width=2200, height=2000, res=400)
    par(mfcol=c(1,1))     
    p <- upset(upSet.df, sets = colnames(upSet.df)[-1], cutoff = NULL,  nintersects = 7,empty.intersections = T)
    print(p)
    dev.off()
}







































































p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_2_1", idents.1 = "2", idents.2 = NULL, bioRep = "name",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 25, outDir = paste0("./output/", outName, "/"), title = "mast_VS_other", idents.1_NAME = "mast", idents.2_NAME = "other", returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 24
                    )

seu.gut.t <- MuDataSeurat::ReadH5AD("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/Tcell_log_counts02_v2.h5ad")
reference <- seu.gut.t

seu.gut.myeloid <- MuDataSeurat::ReadH5AD("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/myeloid_log_counts02_v2.h5ad")
reference <- seu.gut.myeloid

seu.gut.bcell <- MuDataSeurat::ReadH5AD("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/Bcell_log_counts02_v2.h5ad")
reference <- seu.gut.bcell

seu.gut.duod <- MuDataSeurat::ReadH5AD("/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/epi_log_counts02_v2.h5ad")
reference <- seu.gut.duod

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
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_umap_Predicted_tcell.png", sep = ""), width = 10, height = 7)


seu.obj.ros <- readRDS(file = "./output/s3/230828_duod_h3_ros_res0.4_dims50_dist0.4_neigh30_S3.rds")

seu.obj.Hsub <- subset(seu.obj,
                  subset = 
                  cellSource ==  "Healthy"
                 ) 


colnames(seu.obj.ros@assays$RNA@data) <- ifelse(grepl("_1",rownames(seu.obj.ros@meta.data)),paste0(substr(rownames(seu.obj.ros@meta.data), 1, nchar(rownames(seu.obj.ros@meta.data))-2),"_7") ,rownames(seu.obj.ros@meta.data))
colnames(seu.obj.ros@assays$RNA@data) <- ifelse(grepl("_2",rownames(seu.obj.ros@meta.data)),paste0(substr(rownames(seu.obj.ros@meta.data), 1, nchar(rownames(seu.obj.ros@meta.data))-2),"_8") ,rownames(seu.obj.ros@meta.data))
colnames(seu.obj.ros@assays$RNA@data) <- ifelse(grepl("_3",rownames(seu.obj.ros@meta.data)),paste0(substr(rownames(seu.obj.ros@meta.data), 1, nchar(rownames(seu.obj.ros@meta.data))-2),"_9") ,rownames(seu.obj.ros@meta.data))


seu.obj.Hsub <- AddMetaData(seu.obj.Hsub, metadata = seu.obj.ros@assays$RNA@data[rownames(seu.obj.ros@assays$RNA@data) == "TRDC",], col.name = "TRDCexp")

features = c("TRDCexp")
titles <- "TRDC (constant region of TCR delta chain)"
p <- prettyFeats(seu.obj = seu.obj.Hsub, nrow = 1, ncol = 1, features = features,title.size  =14, color = "black", order = F, titles = titles) 
ggsave(paste("./output/", outName,"/", subname, "/", outName, "_supp_TRDC_UMAP.png", sep = ""), width = 7, height = 7)



### Fig supp: plot enrichment scores
ecLists <- read.csv("gut_ecTerms.csv", header = T)

modulez <- split(ecLists$genes, ecLists$cluster)

modulez <- modulez[unname(unlist(lapply(unlist(lapply(modulez, length)), function(x){ifelse(x >= 10, TRUE, FALSE)})))]

names(modulez) <- paste0(names(modulez),"_SIG")

seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")

names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

features <- names(modulez)

fig_supp <- majorDot(seu.obj = seu.obj, groupBy = "majorID",
                     features = features
                    ) + theme(legend.position = "bottom",
                              axis.title.y = element_blank(),
                              plot.margin = margin(7, 7, 0, 100, "pt")) + scale_y_discrete(position = "right") + guides(size = guide_legend(nrow = 2, byrow = F, title = 'Percent\nenriched')) + guides(color = guide_colorbar(title = 'Scaled\nenrichment score')) 
ggsave(paste("./output/", outName, "/", subname, "/", outName, "_gut_modScores.png", sep = ""), width = 12, height = 3, scale=2)


modulez <- list("Mast cell" = c("CPA3","TPSB2","HDC","SLC8A3","SLC18A2","KIT","SLC24A3","GRAP2","IL18R1","GATA2")) ## gene list from Azimuth Kidney
#complete module scoring
seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")

#correct the naming -- credit to: https://github.com/satijalab/seurat/issues/2559
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

#plot the results -- uses a custom function, so you will need to source the customFeunctions.R file. Alt: can also be visulized with FeaturePlot() or DotPlot()
features <- names(modulez)
ecScores <- majorDot(seu.obj = seu.obj, groupBy = "clusterID_2_1", scale = T,
                     features = features
                    ) + theme(axis.title = element_blank(),
                              #axis.ticks = element_blank(),
                              #legend.justification = "left",
                              #plot.margin = margin(7, 21, 7, 7, "pt")
                              legend.direction = "vertical",
                              legend.position = "right"
                             ) + guides(color = guide_colorbar(title = 'Scaled\nenrichment\nscore')) + guides(size = guide_legend(nrow = 3, byrow = F, title = 'Percent\nenriched'))

ggsave(paste("./output/", outName, "/", outName, "_dots_celltypes.png", sep = ""),width = 10,height=6)





colz_df <- read.csv("./colorID_cie3v4.csv")

# colz_df <- colz_df %>% arrange(majorID) %>% mutate(newCol = gg_color_hue(nrow(colz_df)*3)[ c( rep(FALSE, 2), TRUE ) ] )%>% arrange(clusterID_2_1)
# write.csv(colz_df,"./colorID_cie3v4.csv", row.names = F)

labColz <- rep("black", length(levels(seu.obj$clusterID_2_1)))
labColz[25] <- "white"
nudge_x <- rep(0, length(levels(seu.obj$clusterID_2_1)))
nudge_x[28] <- 0.5
nudge_x[19] <- -0.5
nudge_y <- rep(0, length(levels(seu.obj$clusterID_2_1)))
nudge_y[19] <- 0.5

#plot inital cluster umap
pi <- DimPlot(seu.obj, 
        reduction = "umap", 
        group.by = "clusterID_2_1",
              cols=colz_df$newCol,
        pt.size = 0.25,
        label = TRUE,
        label.box = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = labColz, nudge_x = nudge_x, nudge_y = nudge_y) + NoLegend()
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)



#plot inital cluster umap
pi <- DimPlot(seu.obj, 
        reduction = "umap", 
        group.by = "clusterID_2_1",
        pt.size = 0.25,
        label = TRUE,
        label.box = TRUE
 )
p <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = "black") + NoLegend()
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_raw_parent_UMAP.png", sep = ""), width = 7, height = 7)






### Use conocial markers to ID cells
features <- c("PTPRC","CD3G","IL7R", "GZMB",
              "AIF1", "C1QC", "PAX5","JCHAIN",
              "KIT", "ACTA2","COL1A1","SI","MUC13", "TOP2A"
              )


p_feats <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 4, features = features, color = "black", order = T,returnPlots = T,title.size = 16,pt.size = 0.00000001, min.cutoff = "q10") 
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_featPlots.png", sep = ""), width = 12, height = 15)

axes <- ggplot() + labs(x = "UMAP1", y = "UMAP2") + 
theme(axis.line = element_line(colour = "black", 
                               arrow = arrow(angle = 30, length = unit(0.1, "inches"),
                                             ends = "last", type = "closed"),
                              ),
      axis.title = element_text(colour = "black", size = 16),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
     )


patch <- area()
nrow <- 3
ncol <- 6
counter=0
for (i in 1:nrow) {
    if(i < nrow){
        for (x in 3:ncol) {
            counter = counter+1
            if (counter <= length(p_feats)) {
                patch <- append(patch, area(t = i, l = x, b = i, r = x))
            }
        }
    }else{
        for (x in 1:ncol) {
            counter = counter+1
            if (counter <= length(p_feats)) {
                patch <- append(patch, area(t = i, l = x, b = i, r = x))
            }
        }
    }
}

patch <- append(patch, area(t = nrow, l = 1, b = nrow, r = 1))
patch <- append(patch, area(t = 1, l = 1, b = 2, r = 2))

p_final <- Reduce( `+`, p_feats ) +  axes + p + plot_layout(design = patch) 
ggsave(paste0("./output/", outName, "/", subname,"/",outName, "_fig1a.png"), width = 14, height =7)

### Create violin plots for key feats
features = c("CXCL8", "LYZ", "SELL", 
             "S100A12")

pi <- VlnPlot(
    object = seu.obj,
    pt.size = 0,
    same.y.lims = F,
    group.by = "clusterID_2_1",
    combine = T,
    cols = colz_df$newCol,
    stack = T,
    fill.by = "ident",
    flip = T,
    features = rev(features)
        ) + NoLegend() + theme(axis.ticks.y = element_blank(),
                               axis.text.y = element_blank(),
                               axis.title.x = element_blank())

#plot <- prettyViln(plot = pi, colorData = NULL, nrow = 2, ncol = 4)
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_selectViln.png", sep = ""), width = 8, height =5)


### extra code to highlt major groups
lapply(1:length(levels(seu.obj$majorID)), function(x){
    colz <- rep("grey",length(levels(seu.obj$majorID)))
    colz[x] <- colz_df.sub$newCol[x]
    
    umapHighLight <-DimPlot(seu.obj, 
                            reduction = "umap", 
                            group.by = "majorID",
                            cols = colz,
                            pt.size = 0.5,
                            label = F,
                            label.box = F
                           )

    umapHighLight <- formatUMAP(umapHighLight) + NoLegend()
    namE <- levels(seu.obj$majorID)[x]
    ggsave(paste("./output/", outName, "/", subname,"/",outName,"_", namE, "_UMAP.png", sep = ""), width = 7, height = 7)
})




### Use conocial markers to ID cells
features <- c("PTPRC","CD8A","GZMA", 
              "CD3G","CD4", "DLA-DRA")

namez <- c("PTPRC (CD45)", "CD8A", "GZMA", 
           "CD3G", "CD4", "DLA-DRA (MHCII)")


p <- prettyFeats(seu.obj = seu.obj, nrow = 2, ncol = 3, features = features, color = "black", order = F, titles = namez, noLegend = T) 
ggsave(paste("./output/", outName,"/", subname, "/", outName, "_featPlots_tcell.png", sep = ""), width = 9, height = 6)

### Use conocial markers to ID cells
features <- c("ANPEP","EPCAM","SI")


p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, color = "black", order = F, noLegend = T) 
ggsave(paste("./output/", outName,"/", subname, "/", outName, "_featPlots_epithelial.png", sep = ""), width = 9, height = 3)

### Use conocial markers to ID cells
features <- c("KIT","IGF1",
              "TOP2A","JCHAIN")

p <- prettyFeats(seu.obj = seu.obj, nrow = 2, ncol = 2, features = features, color = "black", order = F, noLegend = T) 
ggsave(paste("./output/", outName,"/", subname, "/", outName, "_featPlots_misc.png", sep = ""), width = 6, height = 6)



### Use conocial markers to ID cells
features <- c("PTPRC","CD3G","CD8A", "CD4",
              "GZMA", "S100A12", "DLA-DRA","FLT3",
              "MS4A1","JCHAIN","TOP2A","GATA3",
               "CD34", "ANPEP","EPCAM", "SI",
              "S100G","RELN", "BCAS1","CHGB"
              )


p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 4, features = features, color = "black", order = F) 
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_featPlots.png", sep = ""), width = 12, height = 15)

features <- c("PTPRC","CD68","IL5RA", "LGR5",
              "BMI1", "BMI1", "HOPX","LRIG1",
              "SELL","IL18BP","TERT","CXCL8",
               "CD34", "ESAM","KIT", "IGF1",
              "CX3CR1","ENOH", "LTF","CHGB"
              )


p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 4, features = features, color = "black", order = F) 
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_featPlots_2.png", sep = ""), width = 12, height = 15)


features <- c("PTPRC","CD3G","CD8A", "CD4",
              "CD68","CXCL8","IL18BP","KIT",
              "ANPEP", "EPCAM", "SI","TOP2A"
              )

p <- prettyFeats(seu.obj = seu.obj, nrow = 3, ncol = 4, features = features, color = "black", order = F, noLegend = T) 
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_featPlots_good.png", sep = ""), width = 12, height = 9)



p <- autoDot(seu.integrated.obj = seu.obj, inFile = "/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/viln/allCells_h3c4/230816_duod_h3c4_NoIntrons_gene_list.csv", groupBy = "clusterID_2_1",
                     MIN_LOGFOLD_CHANGE = 0.5, MIN_PCT_CELLS_EXPR_GENE = 0.1,
                    filterTerm = "ENSCAFG"
                    )
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_autoDot.png", sep = ""), width = 9, height = 15)

# p <- autoDot(seu.integrated.obj = seu.obj, inFile = "/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/viln/allCells_h3c4/230713_majorID_duod_h3c4_NoIntrons_gene_list.csv", groupBy = "majorID",
#                      MIN_LOGFOLD_CHANGE = 0.5, MIN_PCT_CELLS_EXPR_GENE = 0.1,
#                     filterTerm = "ENSCAFG"
#                     )
# ggsave(paste("./output/", outName, "/", subname,"/",outName, "_autoDot_2.png", sep = ""), width = 9, height = 15)


p <- dotPlotBY_TYPE(seu_obj = seu.obj, pwdTOvilnCSVoutput = "/pl/active/dow_lab/dylan/k9_duod_scRNA/analysis/output/viln/allCells_h3c4/230816_duod_h3c4_NoIntrons_gene_list.csv", groupBy = "majorID", namedCols = NULL, database = "clfamiliaris_gene_ensembl", exlcude = "", boxColor = "black"
                          )
ggsave(paste("./output/", outName, "/", subname,"/",outName, "_dotByType.png", sep = ""), width = 15, height = 15, scale = 2)


#plot features
fig1c <- majorDot(seu.obj = seu.obj, groupBy = "majorID",
                  yAxis = c("T cell","Epithelial","Plasma cell","Myeloid","Mast cell","Cycling cell","Fibroblast"),
                  features = c("CD3E", "CCL4", "GZMA","GZMB", "CD8A",
                               "SI", "FABP1", "RBP2", "GUCA2A", "APOA1", 
                               "JCHAIN", "RARRES2", "IGHM", "MS4A1", 
                               "AIF1", "C1QC", "S100A12", "LYZ", "CXCL8", 
                               "KIT", "IGF1", "MS4A2",
                               "TOP2A","MKI67",
                              "ACAT2","TPM2","IGFBP7","SFRP1")
                 ) + theme(axis.title = element_blank(),
                           axis.text = element_text(size = 12))
ggsave(paste("./output/", outName, "/", subname, "/", outName, "_majorDot.png", sep = ""), width =8, height = 4)


### Fig 1e: stats by cie vs healthy
freqy <- freqPlots(seu.obj, method = 2, nrow= 1, groupBy = "majorID", legTitle = "Cell source",refVal = "name2", showPval = F,
              namez = "name2", 
              colz = "colz"
              )

freqy + ggpubr::stat_compare_means(method = "t.test",
                                   method.args = list(var.equal = F),
                                   aes(label = paste0("p = ", ..p.format..)), label.x.npc = "left", label.y.npc = 1,vjust = -1, size = 3)


ggsave(paste("./output/", outName, "/", subname, "/",outName, "_freqPlots.png", sep = ""), width = 8.5, height = 3)


Idents(seu.obj) <- "cellSource"
set.seed(12)
seu.obj.ds <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data$cellSource)))
### Fig xx: stacked bar graph by colorID
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "cellSource", groupBy = "name2", clusters = "majorID_pertyName") +
scale_fill_manual(labels = levels(seu.obj$name2), 
               values = levels(seu.obj$colz)) + theme(axis.title.y = element_blank(),
                                                      axis.title.x = element_text(size = 14),
                                                      axis.text = element_text(size = 12)) 
#+ scale_x_discrete(limits=rev(c("Monocyte","DC","B cell","Neutrophil","CD4 T cell","CD8/NK cell","Eosinophil","Basophil","gd T cell","Cycling T cell","CD34+ unk")),expand = c(0, 0))
ggsave(paste("./output/", outName,"/", subname, "/",outName, "_stackedBar.png", sep = ""), width =7, height = 5)



seu.obj$allCells <- "All cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
createPB(seu.obj = seu.obj, groupBy = "allCells", comp = "cellSource", biologicalRep = "name2", lowFilter = T, dwnSam =F, min.cell = 15,
                     clusters = NULL, outDir = paste0("./output/", outName,"/", subname,"/pseudoBulk/") , grepTerm = "H", grepLabel = c("Healthy","CIE") #improve - fix this so it is more functional
                    )

p <- pseudoDEG(metaPWD = paste0("./output/", outName,"/", subname,"/pseudoBulk/allCells_deg_metaData.csv"), returnDDS = F, 
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("./output/", outName,"/", subname,"/pseudoBulk/"), outName = "allCells", idents.1_NAME = "CIE", idents.2_NAME = "Healthy",
          inDir = paste0("./output/", outName,"/", subname,"/pseudoBulk/"), title = "All cells", fromFile = T, meta = NULL, pbj = NULL, returnVolc = T, paired = F, pairBy = "", 
          minimalOuts = F, saveSigRes = T, filterTerm = "^ENSCAF", addLabs = NULL, mkDir = T
                     )

seu.obj$allCells <- "DGE analysis of all cells"
seu.obj$allCells <- as.factor(seu.obj$allCells)
linDEG(seu.obj = seu.obj, threshold = 1, thresLine = F, groupBy = "allCells", comparision = "cellSource", contrast = c("CIE", "Healthy"),
       outDir = paste0("./output/", outName,"/", subname,"/"), 
       outName = "all_cells", cluster = NULL, labCutoff = 15, noTitle = F,
                   colUp = "red", colDwn = "blue", subtitle = T, returnUpList = F, returnDwnList = F, forceReturn = F, useLineThreshold = F, pValCutoff = 0.01, saveGeneList = T, addLabs = ""
                  )

























































### run cell chat on treated
library(CellChat)
seu.obj <- subset(seu.obj, 
                 subset = cellSource == "CIE")

cnts <- seu.obj@assays$RNA@data
meta <- seu.obj@meta.data
# meta$clusterID_sub <- paste("c_",meta$clusterID_sub, sep="")
cell.use <- rownames(meta)
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "celltype.l3") #  usin gpredicted ct.l2

cellchat@idents <- factor(cellchat@idents, levels = as.character(str_sort(levels(cellchat@idents),numeric = TRUE)))

cellchat@DB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#cellchat@DB <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, "./output/cellchat/cellChatobj_cie__ctl3.rds")


#run cell chat on transfered naive
Idents(seu.obj) <- "cellSource"
seu.obj <- subset(seu.obj, 
                 subset = cellSource == "Healthy")

table(seu.obj$cellSource)

cnts <- seu.obj@assays$RNA@data
meta <- seu.obj@meta.data
# meta$clusterID_sub <- paste("c_",meta$clusterID_sub, sep="")
cell.use <- rownames(meta)
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "celltype.l3") #  using predicted ct.l2

cellchat@idents <- factor(cellchat@idents, levels = as.character(str_sort(levels(cellchat@idents),numeric = TRUE)))

cellchat@DB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#cellchat@DB <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, "./output/cellchat/cellChatobj_healthy__ctl3.rds")



#### do the comparision analysis  --- this analysis may be confounded by differences in filter btwn tx and naive samples
outName <- "cellchat"
subName <- ""

cellchat.norm <- readRDS("./output/cellchat/cellChatobj_healthy__ctl3.rds")
cellchat.cie <- readRDS("./output/cellchat/cellChatobj_cie__ctl3.rds")

object.list <- list(Normal = cellchat.norm, OA = cellchat.cie)
object.list[[1]] <- computeCommunProbPathway(object.list[[1]])
object.list[[2]] <- computeCommunProbPathway(object.list[[2]])
object.list[[1]] <- netAnalysis_computeCentrality(object.list[[1]], slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
object.list[[2]] <- netAnalysis_computeCentrality(object.list[[2]], slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

outfile <- paste0("./output/", outName, "/", subName, "/cellchat_cieVn_int.png")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2)) + 
    theme(title = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 20),
          #axis.title.y = element_blank(),
          axis.text = element_text(size = 16)
    ) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1")) + coord_flip()
#gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
#gg1 + gg2
ggsave(outfile, height=2,width=6)

#pretty worthless
outfile <- paste0("./output/", outName, "/", subName, "/cellchat_cieVn_net.png")
png(file = outfile, width=1000, height=500)
par(mfrow = c(1,2), xpd=TRUE)
gg1 <- netVisual_diffInteraction(cellchat, weight.scale = T)
gg2 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()


#make heatmap -- pretty cool
outfile <- paste0("./output/", outName, "/", subName, "/cellchat_cieVn_heatMAP_weight.png")
png(file = outfile, width=2500, height=2500, res=400)
# gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight",
                         cluster.rows = T,
                         cluster.cols = T)
gg2
dev.off()


#plot information flow
outfile <- paste0("./output/", outName, "/", subName, "/cellchat_cieVn_infoFLOW.png")
png(file = outfile, width=1000, height=500)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()

###
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

outfile <- paste0("./output/", outName, "/", subName, "/cellchat_cieVn_int2D.png")
png(file = outfile, width=1000, height=500)
par(mfrow = c(1,2), xpd=TRUE)
gg[[1]] + gg[[2]]
dev.off()




#extract data and customize plots
gg1.df <- gg[[1]]$data
gg1.df$data_type <- "Healthy"

gg2.df <- gg[[2]]$data
gg2.df$data_type <- "CIE"

gg.df <- rbind(gg1.df,gg2.df)

# colz.df <- as.data.frame(colz)
# colz.df$labels <- rownames(colz.df)

# gg.df <- gg.df %>% mutate(strength = x*y) %>% left_join(colz.df, by = "labels")

#set cell types to label
gg.df <- gg.df %>% group_by(data_type) %>% arrange(desc(Count)) %>% mutate(lab=ifelse(row_number() <= 5, as.character(labels), NA)) %>% ungroup()

#make fig 4b
pis <- lapply(c("Healthy","CIE"),function(z){
    gg.df.sub <- gg.df %>% filter(data_type == z)

    ggplot(data=gg.df.sub, aes(x = x, y = y, size=Count, colour = labels, label=lab)) + 
            ggtitle(z) +
            geom_point() + 
            labs(x = "Outgoing interaction strength", y = "Incoming interaction strength") +
            geom_text_repel(max.overlaps = Inf, size=3, color = "black") + 
            theme_classic() + 
    theme(axis.title = element_text(size= 10),
          axis.text = element_text(size= 8),
          title = element_text(size= 11),
          legend.title=element_text(size=10), 
          legend.text=element_text(size=8)
                 ) + NoLegend()
})

pi2 <- Reduce( `+`, pis ) + plot_layout(ncol = 2, guides = 'collect') & scale_size_continuous(limits = c(min(gg.df$Count), max(gg.df$Count)))# & xlim(0, 12.5) & ylim(0, 12.5)
# pi <- pi1 + pi2 + plot_layout(ncol = 2, widths = c(0.22,0.78)) 
ggsave(pi2, file = paste0("./output/", outName, "/", outName, "_interactionScater.png"), width = 6.25, height = 3)

