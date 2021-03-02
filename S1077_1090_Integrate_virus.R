rm(list = ls())

# Prepare working enviroment
library(easypackages)

library("Seurat")
library("RColorBrewer")
library("viridis")
library("xlsx")
library("ggplot2")
library("dplyr")
library("tidyr")
library("pheatmap")
library("tibble")
library(sctransform)
#options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
getwd()
setwd("E:/projects/HT2020-10760/newreference/merge_S1077_1090_virus")
dir()

## Load the Nor, Pre and Ca data of Pt1
S1077_CO_human <- readRDS("E:/projects/analysis_Seurat/newreference/S1077_CO_virus/S1077_CO_human.rds")
S1077_RE_human <- readRDS("E:/projects/analysis_Seurat/newreference/S1077_RE_virus/S1077_RE_human.rds")
S1090_CO_human <- readRDS("E:/projects/analysis_Seurat/newreference/S1090_CO_virus/S1090_CO_human.rds")
S1090_RE_human <- readRDS("E:/projects/analysis_Seurat/newreference/S1090_RE_virus/S1090_RE_human.rds")
## Set up the condition state for each object
S1077_CO_human$condition <- "S1077_CO_human"
S1077_RE_human$condition <- "S1077_RE_human"
S1090_CO_human$condition <- "S1090_CO_human"
S1090_RE_human$condition <- "S1090_RE_human"

## Perform integration
S1077.1090.anchors <- FindIntegrationAnchors(object.list = 
                                        list(S1077_CO_human, S1077_RE_human, S1090_CO_human, S1090_RE_human), 
                                        dims = 1:50)

S1077.1090.combined <- IntegrateData(anchorset = S1077.1090.anchors, dims = 1:50)
dim(S1077.1090.combined) #[1] 2000 9508:共9508个细胞
#查询该整合后的seurat对象的condition种类，以及每个种类所包含的细胞数目
table(S1077.1090.combined$condition)
#S1077_CO_human S1077_RE_human S1090_CO_human S1090_RE_human 
#2827            338           5912            431

###################################
# S1077.1090.combined <- readRDS("E:/projects/analysis_Seurat/virus/merge_S1077_S1090_virus/S1077.1090.combined_jia.rds")
# S1077.1090.combined <- subset(S1077.1090.combined, downsample = 100)
#################################以下注释掉的操作由于还未重新聚类，所以无意义###
# colourCount = length(levels(Idents(S1077.1090.combined)))
# colourCount #[1] 18
# 
# VlnPlot(S1077.1090.combined, features = c("nFeature_RNA", "nCount_RNA","human.percent.mt", 
#                                      "human.percent.RPS","percent.virus"), pt = 0)
# 
# VlnPlot(S1077.1090.combined, features = "percent.virus",ncol = 1, pt = 0.1,log = T)
# 
# VlnPlot(S1077.1090.combined, features = "percent.virus",ncol = 1, pt = 0.1,split.by = 'condition',log = T)


## Perform an integrated analysis
DefaultAssay(S1077.1090.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
S1077.1090.combined <- ScaleData(S1077.1090.combined, verbose = T)

S1077.1090.combined <- RunPCA(S1077.1090.combined, verbose = T)
ElbowPlot(S1077.1090.combined, n = 50)

# t-SNE and Clustering
S1077.1090.combined <- FindNeighbors(S1077.1090.combined, reduction = "pca", dims = 1:30)
S1077.1090.combined
dim(S1077.1090.combined) #[1] 2000 9508
S1077.1090.combined <- FindClusters(S1077.1090.combined, resolution = 1.0)
#resolution:cluster num
#1.5: 20
#1.2：17
#1.1：16
#1.0：16
#0.9：14
S1077.1090.combined <- RunUMAP(S1077.1090.combined, reduction = "pca", dims = 1:30)
S1077.1090.combined <- RunTSNE(S1077.1090.combined, reduction = "pca", dims = 1:30)

##########重新聚类后检查聚成了多少个类###
colourCount = length(levels(Idents(S1077.1090.combined)))
colourCount #[1] 16 ：也就是说聚成了16个cluster
########检查基因等的分布情况###
VlnPlot(S1077.1090.combined, features = c("nFeature_RNA", "nCount_RNA"),
        ncol = 2, pt = 0.1)
VlnPlot(S1077.1090.combined, features = c("human.percent.RPS", "human.percent.mt"), ncol = 2, pt = 0.1)
#病毒分布情况，后边还有一个病毒在每个cluster中的分布情况
##未分组：不添加log与添加log
VlnPlot(S1077.1090.combined, features = c("percent.virus"),
        ncol = 1, pt = 0.1)
VlnPlot(S1077.1090.combined, features = c("percent.virus"), log = T,
        ncol = 1, pt = 0.1)
##分组:不添加log和添加log
VlnPlot(S1077.1090.combined, features = c("percent.virus"), split.by = 'condition',
        ncol = 1, pt = 0.1)
VlnPlot(S1077.1090.combined, features = c("percent.virus"), split.by = 'condition', log = T,
        ncol = 1, pt = 0.1)
#VlnPlot(S1077.1090.combined, features = c("percent.virus"),split.by = 'condition',log = T, cols = getPalette(colourCount),
#        ncol = 1, pt = 0.1)
##### Visualization
table(S1077.1090.combined$condition)

#不分组展示所有细胞的聚类情况
colourCount = length(levels(Idents(S1077.1090.combined)))
mode(colourCount)
colourCount
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))

##不分组,不聚类(tsne法)
DimPlot(S1077.1090.combined, reduction = "tsne", group.by = "condition", 
        pt = 0.1) + NoAxes()
##聚类，不分组，同时加上标签(tsne法+umap法)
DimPlot(S1077.1090.combined, reduction = "tsne", label = T, label.size = 6, repel = F,
        pt = 0.1, cols = getPalette(colourCount)) 
DimPlot(S1077.1090.combined, reduction = "umap", label = T, label.size = 6, repel = F,
        pt = 0.1, cols = getPalette(colourCount)) 
##既分组，又聚类,同时加上标签，但无标签和无坐标轴(tsne法)
DimPlot(S1077.1090.combined, reduction = "tsne", split.by = "condition", 
        pt = 0.1, cols = getPalette(colourCount)) + NoAxes()

##既分组，又聚类,同时加上标签和坐标轴(tsne法+umap法)
DimPlot(S1077.1090.combined, reduction = "tsne", label = T, label.size = 6, repel = F,
        pt = 0.1, cols = getPalette(colourCount), split.by = "condition") 
DimPlot(S1077.1090.combined, reduction = "umap", label = T, label.size = 6, repel = F,
        pt = 0.1, cols = getPalette(colourCount), split.by = "condition") 
#计算condition的每个分组在每个cluster中的数量，并保存入xlsx文件中
S1077.1090.combined_clusters_numbers <- data.frame(table(Idents(S1077.1090.combined), 
                                                         S1077.1090.combined$condition))

write.xlsx(S1077.1090.combined_clusters_numbers, col.names = TRUE, 
           row.names = TRUE, append = FALSE,
           file = "S1077.1090.combined_clusters_numbers.xlsx")

#保存seurat对象为.rds文件
saveRDS(S1077.1090.combined, file = "S1077.1090.combined.rds")

markers <- FindAllMarkers(S1077.1090.combined, assay = "RNA", 
                          logfc.threshold = 0.5, min.pct = 0.5, only.pos = T)

write.xlsx(markers, col.names = TRUE, 
           row.names = TRUE, append = FALSE,
           file = "S1077.1090.combined_clusters_markers.xlsx")

# save.image("E:/projects/analysis_Seurat/newreference/merge_S1077_1090_virus/S1077_1090_Integrate_virus.RData")
# rm(list = ls())


