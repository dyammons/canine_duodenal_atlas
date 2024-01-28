#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")


### Process the seu.obj for coorelation analysis -- percentage of "live cells"

#load annotated object
seu.obj <- readRDS("./output/s3/canine_duodenum_annotated.rds")
outName <- "allCells"

#add new metadata level for coorelation
meta.df <- seu.obj@meta.data
meta.df <- meta.df %>% mutate(majorClass = ifelse(majorID == "epithelial", "Epithelial", "Immune"))
seu.obj@meta.data <- meta.df

#get the true cell type percentages as determined using scRNA seq
groupByList <- seu.obj$name2
clusterList <- seu.obj$majorClass

cluster_freq.table <- as.data.frame(table(groupByList, clusterList)) %>% melt()
cluster_freq.table <- cluster_freq.table[,-3]
colnames(cluster_freq.table) <- c("Sample", "ClusterID", "Count")
cluster_freq.table <- cluster_freq.table %>% dplyr::group_by(Sample) %>% mutate(pct = round(prop.table(Count),2))


### Process flow data - panel 1 t cells / pvt immune vs non-immune

#load flow data
tcell.df <- read.csv("./tcell_flow.csv")

#fix sloppy naming
colnames(tcell.df) <- c("sample","count", "singlets", "live", "immune", "tcell", "cd4_tcell", "cd8_tcell", "dp_tcell", "dn_tcell")
tcell.df$sample <- c("BRUNO","COOP", "FLT6", "HOLSTE", "KYT6", "YYL7")

#calc percentages
tcell.df <- tcell.df %>% mutate(pct_live = live/singlets,
                                
                                pct_live_immune = immune/live,
                                pct_live_nonImmune = (live - immune)/live,
                                
                                pct_live_tcell = tcell/live,
                                pct_immune_tcell = tcell/immune,
                                
                                pct_live_cd4 = cd4_tcell/live,
                                pct_immune_cd4 = cd4_tcell/immune,
                                pct_tcell_cd4 = cd4_tcell/tcell,
                                
                                pct_live_cd8 = cd8_tcell/live,
                                pct_immune_cd8 = cd8_tcell/immune,
                                pct_tcell_cd8 = cd8_tcell/tcell,
                                
                                pct_live_dn = dn_tcell/live,
                                pct_immune_dn = dn_tcell/immune,
                                pct_tcell_dn = dn_tcell/tcell,
                                
                                cd4_tcell_frac = cd4_tcell/(cd8_tcell+cd4_tcell),
                                cd8_tcell_frac = cd8_tcell/(cd8_tcell+cd4_tcell)
                               )


#extract summary statistics
tcell.summary.df <- tcell.df[ ,-c(2:10)] %>%
  summarise(across(
    .cols = is.numeric, 
    .fns = list(Mean = mean, MEDIAN = median, SD = sd, MIN = min, MAX = max), na.rm = TRUE, 
    .names = "{col}-{fn}"
    )) %>% tidyr::pivot_longer(cols = everything()) %>% mutate(NAME = gsub("\\-.*","",name),
                                                               STAT = gsub(".*-","",name),
                                                               value = round((value*100), 2)) %>% select(-name) %>% pivot_wider(names_from = "NAME") %>% as.data.frame() 

rownames(tcell.summary.df) <- tcell.summary.df$STAT
write.csv(tcell.summary.df[ ,-1], file = "./output/tcell_summary_stats.csv")

#add column with info regarding % parent
tcell.df <- tcell.df %>% melt() %>% mutate(parent = ifelse(grepl("pct_immune", variable), "immune", ifelse(grepl("pct_live", variable), "live", "NA" )))
tcell.df$value <- round(tcell.df$value, 2)

#rename levels
pctLive_tcell.df <- tcell.df[tcell.df$variable == "pct_live_immune" | tcell.df$variable == "pct_live_nonImmune", ]
pctLive_tcell.df$ClusterID <- factor(pctLive_tcell.df$variable)
pctLive_tcell.df$variable <- NULL
levels(pctLive_tcell.df$ClusterID) <- c("Immune", "Epithelial")

#use metadata to rename with matching name to scrna names
conversion.df <- read.csv("./bulk_decoder.csv")
pctLive_tcell.df <- pctLive_tcell.df %>% left_join(conversion.df, by = c("sample" = "bulk_name"))
pctLive_tcell.df$parent <- NULL
pctLive_tcell.df <- pctLive_tcell.df[ ,c(4,3,2)]
colnames(pctLive_tcell.df) <- c("Sample", "ClusterID", "pct")

#get scrna and flow data in one df
plotData <- pctLive_tcell.df %>% left_join(as.data.frame(cluster_freq.table)[ ,c(1,2,4)], by = c("Sample", "ClusterID"), suffix = c("Flow", "Scrna"))

#run lin regression
lmModel <- lm(data = plotData, pctScrna ~ pctFlow)

eq <- substitute(italic(y) == a + b ~italic(x)*","~~italic(R)^2~"="~r2, 
         list(a = format(unname(coef(lmModel)[1]), digits = 2),
              b = format(unname(coef(lmModel)[2]), digits = 2),
             r2 = format(summary(lmModel)$r.squared, digits = 3)))

### Fig supp xx: make the coor plot with flow
ggplot(plotData, aes(x=pctFlow, y=pctScrna)) + 
                stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linetype = "dashed", color = "grey50") +
#                 geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
                geom_point(aes(color = ClusterID)) +
                scale_x_continuous(limits = c(0, 1),expand = c(0, 0)) +
                scale_y_continuous(limits = c(0, 1),expand = c(0, 0)) +
                labs(x = "Percentage by flow cytometry", y = "Percentage by scRNA") +
                guides(color = guide_legend(title = "Cell type", size = 3, override.aes=list(fill=NA))) +
                geom_text(x = .5, y = .05, label = as.character(as.expression(eq)), parse = TRUE) +
                theme(panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA),
                      legend.background = element_rect(fill = "transparent",colour = NA),
                      legend.key = element_rect(fill = "transparent",colour = NA),
                      panel.grid.major = element_line(color = "gray"), 
                      panel.grid.minor = element_line(color = "gray"),
                      #axis.text = element_blank(), 
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
                      ) + scale_colour_manual(values=c("#645A9F", "#FF755F", "#009DA5"))
                
ggsave(paste("./output/", outName, "/", outName, "_FlowCoorPlot.png.png", sep = ""), width = 5.25, height = 4)


### Process the seu.obj for coorelation analysis -- percentage of "live cells" at higher res

#add new metadata level for coorelation
meta.df <- seu.obj@meta.data
meta.df <- meta.df %>% mutate(bigClass = ifelse(majorClass == "Epithelial", "Non-immune",
                                                ifelse(majorID %in% c("tcell","cycling"), "T cell", 
                                                ifelse(majorID %in% c("myeloid","mast"), "Myeloid", "Other immune"))))
seu.obj@meta.data <- meta.df

#get the true cell type percentages as determined using scRNA seq
groupByList <- seu.obj$name2
clusterList <- seu.obj$bigClass

cluster_freq.table <- as.data.frame(table(groupByList, clusterList)) %>% melt()
cluster_freq.table <- cluster_freq.table[,-3]
colnames(cluster_freq.table) <- c("Sample", "ClusterID", "Count")
cluster_freq.table <- cluster_freq.table %>% dplyr::group_by(Sample) %>% mutate(pct = round(prop.table(Count),2))


### Process flow data - panel 2 major cell types / pct cd11b vs t cell

#load flow data
major.df <- read.csv("./major_flow.csv")

#fix naming
colnames(major.df) <- c("sample","count", "singlets", "live", "immune", "cd11b", "exclude", "tcell", "other")
major.df$sample <- c("BRUNO","COOP", "FLT6", "KYT6", "YYL7")
major.df <- select(major.df, -exclude)

major.df <- major.df %>% mutate(pct_live = live/singlets,
                                
                                pct_live_immune = immune/live,
                                pct_live_nonImmune = (live - immune)/live,
                                
                                pct_live_tcell = tcell/live,
                                pct_immune_tcell = tcell/immune,
                                
                                pct_live_cd11b = cd11b/live,
                                pct_immune_cd11b = cd11b/immune,
                                
                                pct_immune_other = (immune-cd11b-tcell)/immune,
                                pct_live_other = (live-(live - immune)-cd11b-tcell)/live

                               )

major.summary.df <- major.df[ ,-c(2:8)] %>%
  summarise(across(
    .cols = is.numeric, 
    .fns = list(Mean = mean, MEDIAN = median, SD = sd, MIN = min, MAX = max), na.rm = TRUE, 
    .names = "{col}-{fn}"
    )) %>% tidyr::pivot_longer(cols = everything()) %>% mutate(NAME = gsub("\\-.*","",name),
                                                               STAT = gsub(".*-","",name),
                                                               value = round((value*100), 2)) %>% select(-name) %>% pivot_wider(names_from = "NAME") %>% as.data.frame() 

rownames(major.summary.df) <- major.summary.df$STAT
major.summary.df[ ,-1]
write.csv(major.summary.df[ ,-1], file = "./output/summary_stats.csv")

major.df <- major.df %>% melt() %>% mutate(parent = ifelse(grepl("pct_immune", variable), "immune", ifelse(grepl("pct_live", variable), "live", "NA" )))
major.df$value <- round(major.df$value, 2)

pctLive_major.df <- major.df[major.df$variable == "pct_live_nonImmune" | major.df$variable == "pct_live_tcell" | major.df$variable == "pct_live_cd11b" | major.df$variable == "pct_live_other", ]
pctLive_major.df$ClusterID <- factor(pctLive_major.df$variable)
pctLive_major.df$variable <- NULL
levels(pctLive_major.df$ClusterID) <- c("Non-immune","T cell", "Myeloid", "Other immune")

conversion.df <- read.csv("./bulk_decoder.csv")

pctLive_major.df <- pctLive_major.df %>% left_join(conversion.df, by = c("sample" = "bulk_name"))

pctLive_major.df$parent <- NULL
pctLive_major.df <- pctLive_major.df[ ,c(4,3,2)]
colnames(pctLive_major.df) <- c("Sample", "ClusterID", "pct")

plotData <- pctLive_major.df %>% left_join(as.data.frame(cluster_freq.table)[ ,c(1,2,4)], by = c("Sample", "ClusterID"), suffix = c("Flow", "Scrna"))
plotData <- plotData %>% mutate(pctFlow = pctFlow*100,
                                pctScrna = pctScrna*100)
#run lin regression
lmModel <- lm(data = plotData, pctScrna ~ pctFlow)

eq <- substitute(italic(y) == a + b ~italic(x)*","~~italic(R)^2~"="~r2, 
         list(a = format(unname(coef(lmModel)[1]), digits = 2),
              b = format(unname(coef(lmModel)[2]), digits = 2),
             r2 = format(summary(lmModel)$r.squared, digits = 2)))

### Fig supp xx: make the coor plot with flow
p <- ggplot(plotData, aes(x=pctFlow, y=pctScrna)) + 
                stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linetype = "dashed", color = "grey50") +
#                 geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
                geom_point(aes(color = ClusterID)) +
                scale_x_continuous(limits = c(0, 100),expand = c(0, 0)) +
                scale_y_continuous(limits = c(0, 100),expand = c(0, 0)) +
                labs(x = "Percentage by flow cytometry", y = "Percentage by scRNA") +
                guides(color = guide_legend(title = "Cell type", size = 3, override.aes=list(fill=NA))) +
                geom_text(x = 58, y = 5, label = as.character(as.expression(eq)), parse = TRUE) +
                theme(panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA),
                      legend.background = element_rect(fill = "transparent",colour = NA),
                      legend.key = element_rect(fill = "transparent",colour = NA),
                      panel.grid.major = element_line(color = "gray"), 
                      panel.grid.minor = element_line(color = "gray"),
                      #axis.text = element_blank(), 
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
                      ) + scale_colour_manual(values=c("#C89504", "#FF89B3", "#00ABFD", "#9B8EFF"))
                
ggsave(paste("./output/", outName, "/", outName, "_FlowCoorPlot_fig.png", sep = ""), width = 5.25, height = 3.5)


### Process the seu.obj for coorelation analysis -- percentage of "immune cells"

#subset on immune cells
seu.obj <- subset(seu.obj, subset = majorClass == "Immune")

#add new metadata level for coorelation
meta.df <- seu.obj@meta.data
meta.df <- meta.df %>% mutate(immuneClass = ifelse(majorID %in% c("tcell","cycling"), "T cell", ifelse(majorID %in% c("myeloid","mast"), "Myeloid", "Other immune")))
seu.obj@meta.data <- meta.df

#get the true cell type percentages as determined using scRNA seq
groupByList <- seu.obj$name2
clusterList <- seu.obj$immuneClass

cluster_freq.table <- as.data.frame(table(groupByList, clusterList)) %>% melt()
cluster_freq.table <- cluster_freq.table[,-3]
colnames(cluster_freq.table) <- c("Sample", "ClusterID", "Count")
cluster_freq.table <- cluster_freq.table %>% dplyr::group_by(Sample) %>% mutate(pct = round(prop.table(Count),2))


### Process flow data - panel 2 major cell types / pct cd11b vs t cell

#load flow data
major.df <- read.csv("./major_flow.csv")

#fix naming
colnames(major.df) <- c("sample","count", "singlets", "live", "immune", "cd11b", "exclude", "tcell", "other")
major.df$sample <- c("BRUNO","COOP", "FLT6", "KYT6", "YYL7")
major.df <- select(major.df, -exclude)

major.df <- major.df %>% mutate(pct_live = live/singlets,
                                
                                pct_live_immune = immune/live,
                                pct_live_nonImmune = (live - immune)/live,
                                
                                pct_live_tcell = tcell/live,
                                pct_immune_tcell = tcell/immune,
                                
                                pct_live_cd11b = cd11b/live,
                                pct_immune_cd11b = cd11b/immune,
                                
                                pct_immune_other = (immune-cd11b-tcell)/immune
                               )

major.summary.df <- major.df[ ,-c(2:8)] %>%
  summarise(across(
    .cols = is.numeric, 
    .fns = list(Mean = mean, MEDIAN = median, SD = sd, MIN = min, MAX = max), na.rm = TRUE, 
    .names = "{col}-{fn}"
    )) %>% tidyr::pivot_longer(cols = everything()) %>% mutate(NAME = gsub("\\-.*","",name),
                                                               STAT = gsub(".*-","",name),
                                                               value = round((value*100), 2)) %>% select(-name) %>% pivot_wider(names_from = "NAME") %>% as.data.frame() 

rownames(major.summary.df) <- major.summary.df$STAT
major.summary.df[ ,-1]
write.csv(major.summary.df[ ,-1], file = "./output/myeloid_summary_stats.csv")

major.df <- major.df %>% melt() %>% mutate(parent = ifelse(grepl("pct_immune", variable), "immune", ifelse(grepl("pct_live", variable), "live", "NA" )))
major.df$value <- round(major.df$value, 2)

pctLive_major.df <- major.df[major.df$variable == "pct_immune_tcell" | major.df$variable == "pct_immune_cd11b" | major.df$variable == "pct_immune_other", ]
pctLive_major.df$ClusterID <- factor(pctLive_major.df$variable)
pctLive_major.df$variable <- NULL
levels(pctLive_major.df$ClusterID) <- c("T cell", "Myeloid", "Other immune")

conversion.df <- read.csv("./bulk_decoder.csv")

pctLive_major.df <- pctLive_major.df %>% left_join(conversion.df, by = c("sample" = "bulk_name"))

pctLive_major.df$parent <- NULL
pctLive_major.df <- pctLive_major.df[ ,c(4,3,2)]
colnames(pctLive_major.df) <- c("Sample", "ClusterID", "pct")

plotData <- pctLive_major.df %>% left_join(as.data.frame(cluster_freq.table)[ ,c(1,2,4)], by = c("Sample", "ClusterID"), suffix = c("Flow", "Scrna"))

#run lin regression
lmModel <- lm(data = plotData, pctScrna ~ pctFlow)

eq <- substitute(italic(y) == a + b ~italic(x)*","~~italic(R)^2~"="~r2, 
         list(a = format(unname(coef(lmModel)[1]), digits = 2),
              b = format(unname(coef(lmModel)[2]), digits = 2),
             r2 = format(summary(lmModel)$r.squared, digits = 3)))

### Fig supp xx: make the coor plot with flow
ggplot(plotData, aes(x=pctFlow, y=pctScrna)) + 
                stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linetype = "dashed", color = "grey50") +
#                 geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
                geom_point(aes(color = ClusterID)) +
                scale_x_continuous(limits = c(0, 1),expand = c(0, 0)) +
                scale_y_continuous(limits = c(0, 1),expand = c(0, 0)) +
                labs(x = "Percentage by flow cytometry", y = "Percentage by scRNA") +
                guides(color = guide_legend(title = "Cell type", size = 3, override.aes=list(fill=NA))) +
                geom_text(x = .5, y = .05, label = as.character(as.expression(eq)), parse = TRUE) +
                theme(panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA),
                      legend.background = element_rect(fill = "transparent",colour = NA),
                      legend.key = element_rect(fill = "transparent",colour = NA),
                      panel.grid.major = element_line(color = "gray"), 
                      panel.grid.minor = element_line(color = "gray"),
                      #axis.text = element_blank(), 
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
                      ) + scale_colour_manual(values=c("#645A9F", "#FF755F", "#009DA5"))
                
ggsave(paste("./output/", outName, "/", outName, "_FlowCoorPlot_immune.png.png", sep = ""), width = 5.25, height = 4)




#subset on t cells
seu.obj <- subset(seu.obj, subset = majorID %in% c("tcell", "cycling"))

#add new metadata level for coorelation
meta.df <- seu.obj@meta.data
# meta.df <- meta.df %>% mutate(immuneClass = ifelse(celltype.l3 %in% c("Tnaive","Treg"), "Tinf", ifelse(celltype.l3 %in% c("Trm","gdT_1","gdT_2"), "Trm", "T other")))
meta.df <- meta.df %>% mutate(immuneClass = ifelse(celltype.l3 %in% c("Tnaive","Treg", "CD8mem", "NK_T"), "Non-resident", ifelse(celltype.l3 %in% c("Trm","gdT_1","gdT_2"), "Tissue resident", "T other (ILC2, NK, IFN-T)")))
seu.obj@meta.data <- meta.df

#get the true cell type percentages as determined using scRNA seq
groupByList <- seu.obj$name2
clusterList <- seu.obj$immuneClass

cluster_freq.table <- as.data.frame(table(groupByList, clusterList)) %>% melt()
cluster_freq.table <- cluster_freq.table[,-3]
colnames(cluster_freq.table) <- c("Sample", "ClusterID", "Count")
cluster_freq.table <- cluster_freq.table %>% dplyr::group_by(Sample) %>% mutate(pct = round(prop.table(Count),2))


### Process flow data - panel 1 t cells / pvt immune vs non-immune

#load flow data
tsub.df <- read.csv("./t_subsets_flow.csv")

#fix naming
colnames(tsub.df) <- c("sample","count", "singlets", "live", "immune", "cd11b_neg", "CD4inf", "CD8inf", "Trm")
tsub.df$sample <- c("BRUNO","COOP", "FLT6", "KYT6", "YYL7")

tsub.df <- tsub.df %>% mutate(pct_live = live/singlets,
                                
                                pct_live_immune = immune/live,
                                pct_live_nonImmune = (live - immune)/live,

                                pct_live_susT = cd11b_neg/live,
                                pct_immune_susT = cd11b_neg/immune,

                                pct_live_CD4inf = CD4inf/live,
                                pct_immune_CD4inf = CD4inf/immune,
                                pct_susT_CD4inf = CD4inf/cd11b_neg,

                                pct_live_CD8inf = CD8inf/live,
                                pct_immune_CD8inf = CD8inf/immune,
                                pct_susT_CD8inf = CD8inf/cd11b_neg,
                              
                                pct_live_Trm = Trm/live,
                                pct_immune_Trm = Trm/immune,
                                pct_susT_Trm = Trm/cd11b_neg,
                              
                                pct_live_Tinf = (CD4inf + CD8inf) / live,
                                pct_immune_Tinf = (CD4inf + CD8inf) / immune,
                                pct_susT_Tinf = (CD4inf + CD8inf) / cd11b_neg,
                              
                                pct_susT_other = (cd11b_neg-CD4inf-CD8inf-Trm)/immune
                               )

tsub.summary.df <- tsub.df[ ,-c(2:9)] %>%
  summarise(across(
    .cols = is.numeric, 
    .fns = list(Mean = mean, MEDIAN = median, SD = sd, MIN = min, MAX = max), na.rm = TRUE, 
    .names = "{col}-{fn}"
    )) %>% tidyr::pivot_longer(cols = everything()) %>% mutate(NAME = gsub("\\-.*","",name),
                                                               STAT = gsub(".*-","",name),
                                                               value = round((value*100), 2)) %>% select(-name) %>% pivot_wider(names_from = "NAME") %>% as.data.frame() 

rownames(tsub.summary.df) <- tsub.summary.df$STAT
tsub.summary.df[ ,-1]
write.csv(tsub.summary.df[ ,-1], file = "./output/tsub_summary_stats.csv")

tsub.df <- tsub.df %>% melt() %>% mutate(parent = ifelse(grepl("pct_immune", variable), "immune", ifelse(grepl("pct_susT", variable), "tcell", "NA" )))
tsub.df$value <- round(tsub.df$value, 2)

pctT.tsub.df <- tsub.df[tsub.df$variable == "pct_susT_CD4inf" | tsub.df$variable == "pct_susT_CD8inf" | tsub.df$variable == "pct_susT_Trm" | tsub.df$variable == "pct_susT_Trm" | tsub.df$variable == "pct_susT_Tinf" | tsub.df$variable == "pct_susT_other", ]
pctT.tsub.df$ClusterID <- factor(pctT.tsub.df$variable)
pctT.tsub.df$variable <- NULL
levels(pctT.tsub.df$ClusterID) <- c("CD4inf", "CD8inf", "Tissue resident", "Non-resident", "T other (ILC2, NK, IFN-T)")

conversion.df <- read.csv("./bulk_decoder.csv")

pctT.tsub.df <- pctT.tsub.df %>% left_join(conversion.df, by = c("sample" = "bulk_name"))

pctT.tsub.df$parent <- NULL
pctT.tsub.df <- pctT.tsub.df[ ,c(4,3,2)]
colnames(pctT.tsub.df) <- c("Sample", "ClusterID", "pct")
pctT.tsub.df <- pctT.tsub.df[pctT.tsub.df$ClusterID == "Tissue resident" | pctT.tsub.df$ClusterID == "Non-resident" | pctT.tsub.df$ClusterID == "T other (ILC2, NK, IFN-T)", ]

plotData <- pctT.tsub.df %>% left_join(as.data.frame(cluster_freq.table)[ ,c(1,2,4)], by = c("Sample", "ClusterID"), suffix = c("Flow", "Scrna"))

#run lin regression
lmModel <- lm(data = plotData, pctScrna ~ pctFlow)

eq <- substitute(italic(y) == a + b ~italic(x)*","~~italic(R)^2~"="~r2, 
         list(a = format(unname(coef(lmModel)[1]), digits = 2),
              b = format(unname(coef(lmModel)[2]), digits = 2),
             r2 = format(summary(lmModel)$r.squared, digits = 3)))

### Fig supp xx: make the coor plot with flow
ggplot(plotData, aes(x=pctFlow, y=pctScrna)) + 
                stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linetype = "dashed", color = "grey50") +
#                 geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
                geom_point(aes(color = ClusterID)) +
                scale_x_continuous(limits = c(0, 1),expand = c(0, 0)) +
                scale_y_continuous(limits = c(0, 1),expand = c(0, 0)) +
                labs(x = "Percentage by flow cytometry", y = "Percentage by scRNA") +
                guides(color = guide_legend(title = "Cell type", size = 3, override.aes=list(fill=NA))) +
                geom_text(x = .6, y = .05, label = as.character(as.expression(eq)), parse = TRUE) +
                theme(panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA),
                      legend.background = element_rect(fill = "transparent",colour = NA),
                      legend.key = element_rect(fill = "transparent",colour = NA),
                      panel.grid.major = element_line(color = "gray"), 
                      panel.grid.minor = element_line(color = "gray"),
                      #axis.text = element_blank(), 
                      axis.ticks = element_blank(),
                      axis.title = element_text(size = 14),
                      legend.title = element_text(size = 12),
                      legend.text = element_text(size = 10),
                      #plot.title = element_text(hjust = 0.02),
                      #plot.title.position = "plot",
                      plot.title = element_blank(),
                      title = element_text(size= 16),
                      axis.line = element_blank(),
                      legend.position = c(0.32, 0.82),
                      plot.margin = margin(7, 10, 7, 7, "pt"),
                      panel.border = element_rect(color = "black",
                                                  fill = NA,
                                                  size = 2)
                      ) + scale_colour_manual(values=c("red", "blue", "grey40"))
                
ggsave(paste("./output/", outName, "/", outName, "_FlowCoorPlot_Tsub.png.png", sep = ""), width = 4, height = 4)


### Look at Trm/Tnaive ratio in scrna data
cluster_freq.table <- cluster_freq.table[cluster_freq.table$ClusterID != "T other (ILC2, NK, IFN-T)", ]
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
    ggpubr::stat_compare_means(method = "t.test", method.args = list(var.equal = F), aes(label = paste0("p = ", ..p.format..)), label.x.npc = "left", label.y.npc = 1, vjust = -1, size = 3) + 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    theme(panel.grid.major = element_line(color = "grey", size = 0.25),
          #legend.position = "none",
          text = element_text(size = 12) 
          ) + scale_color_manual(labels = unique(seu.obj$name2), values = levels(seu.obj$colz)) + NoLegend()
ggsave(paste("./output/", outName, "/", outName, "_RatioBoxplot.png", sep = ""), width = 2.5, height = 2.5)

