#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")

### Analysis note: 
# This script loads in the fully annotated processed Seurat object (./output/s3/canine_duodenum_annotated.rds)
# then loads the external cell counts determined using flow cytometry (tcell_flow.csv and major_flow.csv) 
# to complete correlation analysis between gold standard (flow cytometry) and scRNA. 
# The script generates Supplemental figure 2e

############################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   Run correlation analysis   ######## <<<<<<<<<<<<<<
############################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

### Process the seu.obj for coorelation analysis -- percentage of "live cells"

#load annotated object
seu.obj <- readRDS("./output/s3/canine_duodenum_annotated.rds")
outName <- "allCells"

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


### Process flow data - panel 2 major cell types

#load flow data
major.df <- read.csv("./major_flow.csv")
colnames(major.df) <- c("sample","count", "singlets", "live", "immune", "cd11b", "exclude", "tcell", "other")
major.df$sample <- c("BRUNO","COOP", "FLT6", "KYT6", "YYL7")
major.df <- select(major.df, -exclude)

#calc percentages
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

#generate csv with summary stats on the percentages
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
write.csv(major.summary.df[ ,-1], file = "./output/summary_stats_flow_pct.csv")

#clean data for plotting
major.df <- major.df %>% melt() %>% mutate(parent = ifelse(grepl("pct_immune", variable), "immune", ifelse(grepl("pct_live", variable), "live", "NA" )))
major.df$value <- round(major.df$value, 2)
pctLive_major.df <- major.df[major.df$variable == "pct_live_nonImmune" | major.df$variable == "pct_live_tcell" | major.df$variable == "pct_live_cd11b" | major.df$variable == "pct_live_other", ]
pctLive_major.df$ClusterID <- factor(pctLive_major.df$variable)
pctLive_major.df$variable <- NULL
levels(pctLive_major.df$ClusterID) <- c("Non-immune","T cell", "Myeloid", "Other immune")

#join the datasets
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

### Fig supp 2e - make the coor plot with flow
p <- ggplot(plotData, aes(x=pctFlow, y=pctScrna)) + 
                stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linetype = "dashed", color = "grey50") +
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
                      axis.ticks = element_blank(),
                      axis.title = element_text(size= 14),
                      plot.title = element_blank(),
                      title = element_text(size= 16),
                      axis.line = element_blank(),
                      panel.border = element_rect(color = "black",
                                                  fill = NA,
                                                  size = 2)
                      ) + scale_colour_manual(values=c("#C89504", "#FF89B3", "#00ABFD", "#9B8EFF"))
ggsave(paste("./output/", outName, "/", outName, "_supp2e.png", sep = ""), width = 5.25, height = 3.5)

############################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   End correlation analysis   ######## <<<<<<<<<<<<<<
############################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<