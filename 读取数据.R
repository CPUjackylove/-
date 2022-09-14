# July/August 2021
# HBC single-cell RNA-seq workshop

# Single-cell RNA-seq analysis - QC

# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

## Read10X函数会为每个细胞自动创建元数据
ctrl_counts <- Read10X(data.dir = "data/ctrl_raw_feature_bc_matrix")

ctrl <- CreateSeuratObject(counts = ctrl_counts,
                           min.features = 100)
##查看元数据
head(ctrl@meta.data)

##### 读取多个样本，可以使用for循环
## paste0函数，用空字符串连接字符串
for (sample in c('ctrl_raw_feature_bc_matrix','stim_raw_feature_bc_matrix')){
      seurat_data <- Read10X(data.dir = paste0("data/", sample))
      seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                       min.features = 100, 
                                       project = sample) ## project参数是元数据里面的样本标识（sample identity）
      assign(sample, seurat_obj)
}

## 分别查看元数据
head(ctrl_raw_feature_bc_matrix@meta.data)
head(stim_raw_feature_bc_matrix@meta.data)

#####下面将两个对象合并到一个单独的Seurat对象中
merged_seurat <- merge(x=ctrl_raw_feature_bc_matrix, y=stim_raw_feature_bc_matrix,
                       add.cell.id = c('ctrl', 'stim')) ## add.cell.id参数在细胞IDs前加前缀

head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

##
View(merged_seurat@meta.data)

# 计算每个细胞的Novelty score，向元数据中添加一列
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / 
                                  log10(merged_seurat$nCount_RNA)

# 计算每个细胞的线粒体比例，向元数据中添加一列
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

##创建元数据dataframe
metadata <- merged_seurat@meta.data
## 向元数据dataframe添加cell IDs，列名为cells
metadata$cells <- rownames(metadata)
# 创建sample列，内容为每个细胞对应的样本名
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells, "^stim_"))] <- "stim"
# 重命名一些列
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
view(metadata)

# 将新的metadata更新到Seurat对象中
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file="data/merged_filtered_seurat.RData")

