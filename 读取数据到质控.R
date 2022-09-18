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
# Read10X的输出是一个稀疏矩阵
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

#####################

# 可视化每个样本的细胞数量
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# 可视化每个细胞的UMI数量（转录本数量）
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# 直方图可视化每个细胞检测到的基因
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# 可视化复杂度
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# 可视化线粒体比例
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

######################## 过滤

###### 细胞水平的过滤
# 根据阈值过滤低质量的细胞 - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))
####### 基因水平的过滤

# 提取counts
# counts是表达矩阵
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# 返回一个逻辑矩阵nonzero，0值为FALSE，非0值为TRUE
nonzero <- counts > 0

# 对每一行（每个基因），统计TRUE的数量（计算行和），大于等于10则返回TURE。
# rowSums()函数返回一维向量
# 整个表达式返回一维逻辑向量
keep_genes <- Matrix::rowSums(nonzero) >= 10

# 只保留在10个及以上细胞中保留的基因
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)



