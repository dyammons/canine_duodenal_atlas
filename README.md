[![DOI](https://zenodo.org/badge/635773561.svg)](https://zenodo.org/doi/10.5281/zenodo.11153534)

# Canine duodenal atlas

This GitHub repository contains all the analysis code used in, "Single-Cell RNA Sequencing Implicates Neutrophil and Epithelial Contributions to the Pathogenesis of Chronic Inflammatory Enteropathy in Dogs."

If you use our raw/processed data, extract data using the UCSC Cell Browser portal, or use portions of our code in your analysis, please cite:
> Manchester, A. C., Ammons, D. T., Dow, S., & Lappin, M. R. Single cell transcriptomic analysis of the canine duodenum in chronic inflammatory enteropathy and health. Frontiers in Immunology, 15, 1397590. doi: 10.3389/fimmu.2024.1397590

## Repository goals: 
- provide a resource to make the data generated from this project accessible
- enable reproducible/transparent data reporting
- provide analysis code to reproduce custom figures

If you have any questions or concerns, please submit an issue, contact the corresponding author(s), and/or contact Dylan Ammons at dylan.ammons @ colostate dot edu.

## File structure:
- [:file\_folder: input](/scripts/input) contains relevant metadata files and instructions for obtaining data associated with this study
- [:file\_folder: scripts](/scripts) contains the analysis code and source file used to complete the data analysis
- [:file\_folder: output](/scripts/output) contains the expected directory structure

## Supplemental data and potential uses:
1. [Browse the data](#1-browse-the-complete-annotated-dataset)
2. [Cell type annotations](#2-cell-type-annotations-with-defining-markers)
3. [Reference Mapping](#3-using-the-data-to-complete-reference-mapping)
4. [Module scoring](#4-module-scoring)

### 1. Browse the complete annotated dataset

The proccessed dataset is avaliable for browsing via the UCSC Cell Browser portal.
Using the portal you can explore feature expression throughout the dataset as well as obtain the transcriptomic signatures of each cell type though an interactive webpage.

Link to the dataset: https://canine-duodenum-cie.cells.ucsc.edu

Link to UCSC Cell Browser documentation: https://cellbrowser.readthedocs.io/en/master/

### 2. Cell type annotations with defining markers

Cell markers lists were curated using the full dataset of 3 healthy and 4 chronic inflammatory enteropathy (CIE) dogs. The top 50 defining features (identified using `FindAllMarkers` for each cell type were considered, with the top 24 features evaluated for specificity using violin plots and preference given to unique features only found in the top 50 of one cell type.

<details open><summary>Cell type gene signatures</summary>
<p>

|Cell type        |                                               |Major markers                                                 |
|-----------------|-----------------------------------------------|--------------------------------------------------------------|
|T cells          |                                               |                                                              |
|                 |Effector CD8 (CD8_eff)                         |GZMA, GZMB, CCL5, CD7, CD2, ITGB7, CD96, CD3E                 |
|                 |Naïve (Tnaive)                                 |VIM, IL7R, S100A8, EEF1A1, CXCR4, LTB, TMSB10                 |
|                 |Tissue resident memory (CD8_TRM)               |CCL4, RGS1, IFNG, NR4A2, BCL2A1, FASLG                        |
|                 |gd T subset 1 (gdT_1)                          |GZMA, TRAT1, GZMB, TNFRSF6B, PTPN22, CCL5                     |
|                 |Natural killer T (NK_T)                        |CRTAM, CTSW, CD160, TCF7, TESC, NME1, REL                     |
|                 |Memory CD8 (CD8_mem)                           |GZMK, FOS, SH2D1A, CCL4, CXCR4, KLF4, DLA-DRA                 |
|                 |gd T subset 2 (gdT_2)                          |ENSCAFG00000030206, CCL5, SPNS3, GZMA, CAPG, IL2RB, TRAT1, CD7|
|                 |Regulatory T (Treg)                            |TNFRSF4, TNFRSF18, S100A5, CTLA4, ICOS, NFKBIA                |
|                 |Natural killer (NK)                            |KRLB1, NFKBID, NCR3, F2RL3, IL12RB2, CTSW                     |
|                 |Interferon-enriched T (T_IFN)                  |ISG15, MX2, IFGGB2, OAS1, IFI44, SAMD9L, IFI44L               |
|                 |type 2 innate lymphoid (ILC2)                  |MMP9, IL17RB, GATA3, PLAC8A, CLINT1, MPP4, ARL4C              |
|                 |CD20+ T (MS4A1_T_cell)                         |MS4A1, CD52, TESC, FCER1G, STAP1, GZMB, IL2RB                 |
|Epithelial cells |                                               |                                                              |
|                 |Enterocyte 1                                   |ND4L, UGT1A6, SI, MYO1A, FAM13A, FAM15A, FAM151A              |
|                 |Enterocyte 2                                   |APOC3, APOA4, APOA1, APOB, TM4SF20, CLCA4, SLC40A1            |
|                 |Enterocyte 3                                   |FABP1, ATP5MC1, MGST3, COX4I1, RBP2, GSTP1, SLC25A5           |
|                 |BEST4+ enterocyte                              |GUCA2A, GUCA2B, CFTR, SYNC, STOM, CA7, BEST4                  |
|                 |Goblet cell                                    |ZG16, SPINK4, CLCA1, BCAS1, LYPD8, PNLIP, SYTL2               |
|                 |Tuft cell                                      |ANXA4, CA2, TRPM5, IRAG2, RYR1, FYB1, SYNJ1, ATP2A3           |
|                 |Interferon-enriched enterocyte (IFN_enterocyte)|ISG15, RSAD2, IFIT1, DDX60, APOC3, RNF213, OAS1               |
|                 |Stromal cell                                   |TPM2, VIM, ADAMDEC1, IGFBP7, COL1A2, MYL9                     |
|                 |Enteroendocrine cell                           |CHGB, ADGRG4, SCG2, AFP, UNC13B, PCSK1N, TPH1                 |
|Myeloid cells    |                                               |                                                              |
|                 |Neutrophil                                     |S100A12, S100A8, SERPINA1, SOD2, NFKBIA, SELL, ISG20          |
|                 |Eosinophil                                     |MMP1, MS4A2, CAT, PGF, CHI3L1, ADAMDEC1, FABP3                |
|                 |Monocyte                                       |MT1E, MT2A, HMOX1, C1QC, C1QA, GSTP1, CTSS                    |
|                 |Macrophage                                     |CCL3, MAFB, C1QA, C1QC, STAB1, CTSS, CCL7, C1QB               |
|                 |cDC1                                           |TMSB10, BATF3, ECRG4, CLEC1B, NAPSA, DNASE1L3                 |
|                 |IL22RA2 DC                                     |IL22RA2, TMBS10, CD52, ECRG4, MAL, S100A5, FBP1, CA2, LSP1    |
|Miscellaneous    |                                               |                                                              |
|                 |Plasma cells                                   |JCHAIN, IGHM, TXNDC5, RARRES2, PRDX4, DERL3                   |
|                 |Cycling T cells                                |TUBA1B, TOP2A, STMN1, CENPF, HMGB2                            |
|                 |Mast cells                                     |KIT, IGF1, CD52, CXCR4, SVIL, F2RL3                           |
|                 |B cells                                        |MS4A1, TNFRSF13C, DLA-DRA, FCRLA, SAMSN1, CD22, PAX5          |


</p>
</details>

### 3. Using the data to complete reference mapping
Reference mapping is useful tool to facilitate the identification of cell types in single cell datasets. The approach described here uses Seurat functions to identify anchors between a query dataset (external/personal data) and the reference datasets generated in this study. The default approach describes how to use the healthy only dataset, but it will also work with the combined dataset if you load that file in as the reference.

Before running the reference mapping code, a Seurat object need to be preprocessed and stored as an object named `seu.obj`.
```r
#set the path to the location in which the reference file is saved
reference <- readRDS(file = "AllCells_duod_annotated.rds")

#prepare the reference
reference[['integrated']] <- as(object = reference[['integrated']] , Class = "SCTAssay")
DefaultAssay(reference) <- "integrated"

#find conserved anchors with query and reference
anchors <- FindTransferAnchors(
    reference = reference,
    query = seu.obj,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims= 1:50
)

#select meta.data slot to use for label transfer -- change refdata value to use alternate labels (i.e., refdata = reference$celltype.l1)
predictions <- TransferData(
    anchorset = anchors,
    refdata = reference$celltype.l3,
    dims = 1:50
)
seu.obj <- AddMetaData(seu.obj, metadata = predictions)

#generate and save the image
pi <- DimPlot(
    seu.obj, 
    reduction = "umap", 
    group.by = "predicted.id",
    pt.size = 0.25,
    label = T,
    label.box = T,
    shuffle = F
)
ggsave("./output/referenceMap.png", width = 7, height = 7)
```

### 4. Module scoring
Module scoring is a supplemental approach that can be applied to single cell datasets with the goal of providing further insights into cell identities. The approach described below uses the Seurat function `AddModuleScore` and the gene lists presented in Table 2 of our associated manuscript. 

The concept of the AddModuleScore() function is similar to GSEA, but also distinct in many ways. Read the [Seurat documentation](https://satijalab.org/seurat/reference/addmodulescore) and/or check out [this webpage](https://www.waltermuskovic.com/2021/04/15/seurat-s-addmodulescore-function/) for more details.

```r
#load in the reference file
# can dowlaned with this command: wget https://raw.githubusercontent.com/dyammons/canine_duodenal_atlas/main/scripts/input/genesig_long.csv
ref.df <- read.csv("genesig_long.csv", header = T)

#organize the data
modulez <- split(ref.df$gene, ref.df$celltype.l2)

#complete module scoring
seu.obj <- AddModuleScore(seu.obj, features = modulez, name = "_score")

#correct the naming -- credit to: https://github.com/satijalab/seurat/issues/2559
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

#plot the results
features <- names(modulez)
ecScores <- DotPlot(
    seu.obj,
    assay = "RNA",
    features = features
)

ggsave("./output/dots_celltypes.png", width = 10, height = 6)
```


