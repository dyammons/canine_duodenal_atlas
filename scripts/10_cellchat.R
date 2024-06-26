#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")
library(CellChat)

### Analysis note: 
# This script loads in annotated object and completes cellchat analysis 

############################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin cell chat analysis   ######## <<<<<<<<<<<<<<
############################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

seu.obj <- readRDS("./output/s3/canine_duodenum_annotated.rds")
outName <- "cellchat"

#extract data for cellchat
seu.obj.backup <- seu.obj

### Run cell chat on CIE samples
seu.obj <- subset(seu.obj, 
                  subset = cellSource == "CIE")

#get 1:1 orthologues
cnts <- seu.obj@assays$RNA@data
cnts <- orthogene::convert_orthologs(gene_df = cnts,
                                     gene_input = "rownames", 
                                     gene_output = "rownames", 
                                     input_species = "dog",
                                     output_species = "human",
                                     non121_strategy = "drop_both_species") 
rownames(cnts) <- unname(rownames(cnts))

meta <- seu.obj@meta.data
cell.use <- rownames(meta)
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "celltype.l3")

cellchat@idents <- factor(cellchat@idents, levels = as.character(str_sort(levels(cellchat@idents),numeric = TRUE)))
cellchat@DB <- CellChatDB.human

#run standard cellchat v1 workflow
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, "./output/cellchat/Apr_28_2024_cellChatobj_cie_ctl3.rds")


### Run cell chat on healthy samples
seu.obj <- seu.obj.backup

Idents(seu.obj) <- "cellSource"
seu.obj <- subset(seu.obj, 
                 subset = cellSource == "Healthy")

#get 1:1 orthologues
cnts <- seu.obj@assays$RNA@data
cnts <- orthogene::convert_orthologs(gene_df = cnts,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "dog",
                                        output_species = "human",
                                        non121_strategy = "drop_both_species") 
rownames(cnts) <- unname(rownames(cnts))

meta <- seu.obj@meta.data
cell.use <- rownames(meta)
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "celltype.l3")

cellchat@idents <- factor(cellchat@idents, levels = as.character(str_sort(levels(cellchat@idents),numeric = TRUE)))
cellchat@DB <- CellChatDB.human

#run standard cellchat v1 workflow
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, "./output/cellchat/Apr_28_2024_cellChatobj_healthy_ctl3.rds")



###############################
### BEGIN CELLCHAT ANALYSIS ###
###############################

#set output specifications
outName <- "cellchat"
subName <- "cie_vs_h"

#restore full Seurat object
seu.obj <- seu.obj.backup
# seu.obj <- readRDS("./output/s3/Dec_19_2023_allCells_2500_res0.7_dims40_dist0.3_neigh30_S3.rds") #alterantively reload
# seu.obj.backup <- seu.obj

#load in processed cellchat data
# cellchat.norm <- readRDS("./output/cellchat/Nov_30_2023_cellChatobj_norm_ctl2.rds")
cellchat.norm <- readRDS("./output/cellchat/Apr_28_2024_cellChatobj_healthy_ctl3.rds")
cellchat.cie <- readRDS("./output/cellchat/Apr_28_2024_cellChatobj_cie_ctl3.rds")

colz.base <- gg_color_hue(length(levels(seu.obj$celltype.l3)))
names(colz.base) <- levels(seu.obj$celltype.l3)

#get overlapping pathways
print(cellchat.cie@netP$pathways[cellchat.cie@netP$pathways %in% cellchat.norm@netP$pathways])

#merge cellchat objects
object.list <- list(Healthy = cellchat.norm, CIE = cellchat.cie)
object.list[[1]] <- computeCommunProbPathway(object.list[[1]])
object.list[[2]] <- computeCommunProbPathway(object.list[[2]])
object.list[[1]] <- netAnalysis_computeCentrality(object.list[[1]], slot.name = "netP")
object.list[[2]] <- netAnalysis_computeCentrality(object.list[[2]], slot.name = "netP") 
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

### Fig supp 4a - bar chart with number of interactions
outfile <- paste0("./output/", outName, "/", subName, "_supp4a.png")
data.df <- compareInteractions(cellchat, show.legend = F)$data

data.df$Dataset <- factor(data.df$dataset, levels = c("CIE", "Healthy"))
p <- ggplot(data.df, aes(x=count, y=Dataset, fill=dataset)) + geom_bar(stat="identity") + theme_void() + 
    theme(title = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 18),
          axis.title.y = element_blank(),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(size = 16, hjust = 1),
          axis.text.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 0)),
          plot.margin = margin(4, 14, 4,4)
    ) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1"))  + scale_y_discrete(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) + NoLegend() + xlab("Number of inferred interactions")
ggsave(outfile, height = 1, width = 6)


### Fig 4a - interactivity scater plots by condition
#create plot using CellChat built in function
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

outfile <- paste0("./output/", outName,"//",outName,"_cellchat_cieVh_int2D.png")
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

colz.df <- as.data.frame(colz.base)
colz.df$labels <- rownames(colz.df)

gg.df <- gg.df %>% mutate(strength = x*y) %>% left_join(colz.df, by = "labels")

#set cell types to label
gg.df <- gg.df %>% group_by(data_type) %>% arrange(desc(strength)) %>% mutate(lab=ifelse(row_number() <= 5, as.character(labels), NA)) %>% ungroup()

#make figure
pis <- lapply(c("Healthy","CIE"),function(z){
    gg.df.sub <- gg.df %>% filter(data_type == z)

    ggplot(data=gg.df.sub, aes(x = x, y = y, size=Count, colour = labels, label=lab)) + 
            ggtitle(z) +
            geom_point(color = gg.df.sub$colz.base) + 
            labs(x = "Outgoing interaction strength", y = "Incoming interaction strength") +
            geom_text_repel(max.overlaps = Inf, size = 3, color = "black", box.padding = 0.5,  min.segment.length = 0.2, max.iter = 100000, seed = 666) + 
            theme_classic() + 
    theme(axis.title.x = element_text(size= 10),
          axis.title.y = element_text(size= 10),
          axis.text = element_text(size=8),
          title = element_text(size= 11),
          legend.title = element_text(size=10), 
          legend.text = element_text(size=8),
          legend.position = "top",
          legend.direction = "horizontal",
          plot.title = element_text(face = "bold", hjust = 0.5)
                 )
})

pis[1][[1]] <- pis[1][[1]] + theme(axis.title.x = element_blank(),
                                   axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank()) + guides(size = guide_legend(nrow = 1))

pi2 <- Reduce( `+`, pis ) + plot_layout(ncol = 1, guides = 'collect') & 
    scale_size_continuous(limits = c(min(gg.df$Count), max(gg.df$Count))) & 
    xlim(0, 3.5) & ylim(0, 3.5) & 
    scale_color_manual(values = unname(colz.base), labels = names(colz.base)) & 
    theme(legend.position = "bottom")
ggsave(pi2, file = paste0("./output/", outName, "/", subName, "_fig4a.png"), width = 3.25, height = 6)


### Fig 4b - create differential interaction heatmap
#create heatmap using CellChat built-in function
gg2 <- netVisual_heatmap(cellchat, measure = "weight",
                         cluster.rows = F,
                         cluster.cols = F) # , color.use = colz.base

#extract the data and customize the heatmap
colz.base <- colz.base[match(rownames(gg2@matrix), names(colz.base))]
ha <- HeatmapAnnotation(celltype = names(colz.base),
                        col = list(celltype = colz.base),
                        simple_anno_size = grid::unit(0.2, "cm"),
                        show_annotation_name = FALSE,
                        show_legend = FALSE
                       )
row_ha <- rowAnnotation(celltype = names(colz.base),
                        col = list(celltype = colz.base),
                        simple_anno_size = grid::unit(0.2, "cm"),
                        show_annotation_name = FALSE,
                        show_legend = FALSE
                       )

outfile <- paste0("./output/", outName, "/", subName, "_fig4b.png")
png(file = outfile, width=2800, height=2500, res=500)

cusHeat <- Heatmap(
    gg2@matrix, 
    na_col = "white", 
    col = gg2@matrix_color_mapping,
    bottom_annotation = ha, 
    left_annotation = row_ha,
    show_column_dend = FALSE, 
    show_row_dend = FALSE, 
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    column_title = "Differential interaction strength",
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot = 60,
    row_title = "Sources (Sender)", 
    row_title_gp = gpar(fontsize = 10),
    row_title_rot = 90,
    cluster_rows=T, 
    cluster_columns=T,
    heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"), title_position = "leftcenter-rot",
                                border = NA,
                                legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
)

cusHeat
dev.off()

### Fig 4c - stacked bargraph of network information flow
#plot information flow using CellChat built-in function
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)

#extract data and plot using custom code
res.df <- gg1$data %>% filter(pvalues < 0.05)

orderz <- res.df[res.df$group == "Healthy", ] %>% left_join(res.df[res.df$group == "CIE", ], by = "name") %>% mutate(pct = contribution.x/contribution.y) %>% pull(name) %>% as.character()

res.df$name <- factor(res.df$name, levels = orderz)
res.df$group <- factor(res.df$group, levels = c("Healthy", "CIE"))
p <- ggplot(res.df, aes(x = contribution, y = name, fill = factor(group))) +
            geom_bar(stat = "identity", position = "fill", width = 1, colour="white") +
            theme_classic() + geom_vline(xintercept = 0.5, linetype="dashed", color = "grey50", size=1) +
            theme(title = element_text(size = 14),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 16),
                  legend.position = "top",
                  legend.direction = "horizontal",
                  plot.title = element_blank(),
                  axis.line = element_line(colour = "black"),
                  legend.key.size = unit(1,"line"),
                  axis.title = element_text(size = 18),
                  axis.text = element_text(size = 16),
                  plot.margin = margin(t = 0, r = 21, b = 0, l = 0, unit = "pt")
            ) +
            scale_y_discrete(expand = c(0, 0)) +
            scale_x_continuous(expand = c(0,0)) + 
            ylab(label = "Interaction network") +
            xlab(label = "Relative information flow") + 
            guides(fill = guide_legend(nrow = 1)) + scale_fill_manual(values=c("mediumseagreen","mediumpurple1"))

ggsave(p, file = paste0("./output/", outName, "/", subName, "_fig4c.png"), height = 8, width = 7)



#fix annotation error
groupz.df <- table(seu.obj$majorID, seu.obj$celltype.l3) %>% melt() %>% filter(value > 0)
groupzNames <- groupz.df$Var.1
names(groupzNames) <- groupz.df$Var.2

cellchat <- cellchat.cie
pathwayz <- c("VEGF","IL1","NOTCH", "ADGRE5")
pathways <- pathwayz
subName <- "cie_vs_h"

cellchat <- cellchat.norm
pathwayz <- c("CD40","LT", "PD-L1", "IL2", "GAS", "CD96", "FASLG")
pathways <- pathwayz
subName <- "sf"

pis <- lapply(pathwayz, function(pathway){
    
    #extract required plotting data
    lrData <- as.data.frame(cellchat@LR)
    net <- cellchat@net
    prob <- net$prob
    prob <- prob[,,rownames(lrData[lrData$LRsig.pathway_name == pathway,])]
    prob.sum <- apply(prob, c(1,2), sum)

    #identify which cell types are active in pathway
    idx1 <- which(Matrix::rowSums(prob.sum) == 0)
    idx2 <- which(Matrix::colSums(prob.sum) == 0)
    idx <- intersect(idx1, idx2)
    net <- prob.sum[-idx, ]
    if(is.matrix(net)){
        net <- net[, -idx]
        cellTypeOFinterest <- rownames(net)
    }else{
        net <- net[-idx]
        cellTypeOFinterest <- names(net)
    }
    
    #conditional to fix problem that occurs if all cell types are active in the plot
    if(is.null(cellTypeOFinterest)){
        cellTypeOFinterest <- levels(cellchat@idents)
    }

    #grey out cell types not involved
    colz2 <- colz.base
    colz2[names(colz2)[!names(colz2) %in% cellTypeOFinterest]] <- "grey"

    #save the plot
    outfile <- paste("./output/", outName, "/", subName, "_", pathway ,"_cell_cell_3.png", sep = "")
    png(file = outfile, width=2500, height=2500, res=400)
    par(mfrow=c(1,1), mar = c(0,0,0,0))
    gg7 <- netVisual_aggregate_mod(cellchat, layout = "chord", signaling = pathway, color.use = colz2, remove.isolate = F, big.gap = 5, group = groupzNames, title.name = NULL)
    dev.off()

    #get the active features in the pathway
    genez <- lapply(pathways, function(x){extractEnrichedLR(cellchat, signaling = x, geneLR.return = TRUE, enriched.only = T)[["geneLR"]]})
    names(genez) <- pathways

    Idents(seu.obj) <- "celltype.l3"
    seu.obj.sub <- subset(seu.obj, idents = cellTypeOFinterest)
    
    #plot expression using Seurat function
    pi <- VlnPlot(
        object = seu.obj.sub,
        pt.size = 0,
        same.y.lims = T,
        flip = T,
        group.by = "celltype.l3",
        fill.by = "ident",
        split.by = "cellSource",
        stack = TRUE,
        combine = FALSE,
        features = unlist(genez[pathway])
    ) + NoLegend() + theme(axis.title.x = element_blank(),
                           axis.title.y.right = element_blank())
    
    ggsave(paste("./output/", outName, "/", subName, "_", pathway ,"_viln.png", sep = ""), height = 7, width = 7)
    
    return(pi)
    
}) 

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end cell chat analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<