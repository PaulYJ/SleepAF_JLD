## Milo
## Milo is available from Bioconductor (preferred stable installation)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c("miloR", "SingleCellExperiment", "scater", "scran", "dplyr", "patchwork")
BiocManager::install(packages)

## Install development version
devtools::install_github("MarioniLab/miloR", ref="devel") 
install.packages("ggplot2")
install.packages("Seurat")
install.packages("igraph")


library(readr)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)
library(ggplot2)

dataset <- read_csv("~/Documents/MATLAB/afib_predict.csv")
dataset <- t(dataset)
data_matrix <- as.matrix(dataset)

variables <- c('slpprdp', 'slpeffp', 'timerem', 'timest1', 'timest2', 'timest34',
               'timeremp', 'timest1p', 'timest2p', 'times34p',
               'hslptawp', 'waso', 'ai_all', 'ai_nrem',
               'ai_rem', 'slplatp', 'remlaiip')

# 提取表达矩阵并删除变量名称中的下划线
expression_matrix <- dataset[variables, ]
expression_matrix <- apply(expression_matrix, 2, as.numeric) # 确保数据为数值型

# 对数据进行标准化（这里使用 R 基本包中的 scale 函数）*注意是对每一行进行标准化
scaled_data <- t(scale(t(expression_matrix)))
print(head(scaled_data[1:4,1:4]))

# 创建 meta.features DataFrame，应包含特征名
feature_names <- rownames(expression_matrix)  # 或你特征的正确名字列表
meta_features <- data.frame(row.names = feature_names)

# 创建一个新的 Seurat 对象
seurat_object <- CreateSeuratObject(counts = expression_matrix)
# 将 scaled_data 添加为一个新层
seurat_object <- SetAssayData(object = seurat_object, slot = "scale.data", new.data = scaled_data, assay = "RNA")
# 确保新层已经被成功添加
print(seurat_object[["RNA"]]@layers)

#添加umap结果
umap_embeddings <- read_csv("~/Documents/MATLAB/umap_embeddings.csv", col_names = c("umap_1", "umap_2"))
umap_matrix <- as.matrix(umap_embeddings)
rownames(umap_matrix) <- paste("Cell", 1:nrow(umap_matrix), sep = "_")
seurat_object[["umap"]] <- CreateDimReducObject(embeddings = umap_matrix, key = "umap_", assay = "RNA")

# 提取元数据（假设剩余列是元数据）并添加到 Seurat 对象
tdataset = as.data.frame(t(dataset))
meta_data = as.data.frame(tdataset)
colnames(expression_matrix) = c(1:7841)

rownames(meta_data) <- colnames(expression_matrix)  # 确保元数据的行名与细胞名称匹配
seurat_object <- AddMetaData(seurat_object, metadata = meta_data)


# 运行 PCA 并调整 npcs 参数
# 添加标准化后的数据为新的 scale.data 层
npcs <- 15  # 你可以根据数据集的大小调整这个值
seurat_object <- RunPCA(seurat_object, npcs = npcs, assay = "RNA", slot = "scale.data")
pca_embeddings <- Embeddings(seurat_object, reduction = "pca")
# 查看每个主成分的解释方差
ElbowPlot(seurat_object, ndims = npcs)
# 可视化前几个主成分的基因负载
VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca")


# 建立milo对象
psg_sce <- as.SingleCellExperiment(seurat_object)
psg_milo <- Milo(psg_sce)
reducedDim(psg_milo, "UMAP") <- reducedDim(psg_sce, "UMAP")

## Constructing kNN graph with k:10
psg_milo <- buildGraph(psg_milo, k = 15, d = 9)

## 用makeNhoods函数定义领域
# prop: 初始随机抽样的细胞比例，一般0.1-0.2.
# k: 用于KNN细化的邻域大小 一般和KNN图创建保持一致
# d: 用于KNN细化的降维数 一般和KNN图创建保持一致
# refined: 是否使用细化算法进行采样 默认推荐使用

set.seed(9) # 设置随机种子确保数字可重复
psg_milo <- makeNhoods(psg_milo, prop = 0.2, k = 15, d = 9, refined = TRUE) 

#可视化邻域大小的分布，通过调整前面k和prop的参数获得合适大小
plotNhoodSizeHist(psg_milo)

# 假设 psg_milo 已经有一个包含年龄的列名为 'age' （创建agegroup）
psg_milo$age <- as.numeric(psg_milo$age)
psg_milo$agegroup <- cut(psg_milo$age,
                         breaks = c(-Inf,55,60,65,70,75,80, Inf),
                         labels = c("<55","55-60", "60-65","65-70", "70-75","75-80", ">80"),
                         right = FALSE)  # right = FALSE 表示左闭右开

#给milo对象添加condition
subtype_uni <- "afib" # 可以将此变量名更改为您希望的任何名称
psg_milo$condition <- ifelse(psg_milo@colData@listData[["afibPredict"]] == 1, "R", "NR")
psg_milo$group <- paste(psg_milo$condition, psg_milo$agegroup, sep = "_")
#计算每个领域所含的细胞数
colnames(colData(psg_milo))
psg_milo <- countCells(psg_milo, meta.data = data.frame(colData(psg_milo)), samples="group")
head(nhoodCounts(psg_milo))
sum(nhoodCounts(psg_milo))

#准备进行差异丰度分析的实验设计矩阵
psg_design <- data.frame(colData(psg_milo))[,c("group", "condition")]
psg_design <- distinct(psg_design)
rownames(psg_design) <- psg_design$group
## 根据邻域计数矩阵的列名，重排行名
# psg_design <- psg_design[match(colnames(nhoodCounts(psg_milo)), rownames(psg_design)), ]

# 差异丰度分析统计检验

da_results <- testNhoods(psg_milo, design = ~ condition, design.df = psg_design)

table(da_results$SpatialFDR < 0.1)

# 查看统计量
da_results %>%
  arrange(SpatialFDR) %>%
  head() 

#检查 DA 测试结果
# 首先检查未校正 P 值的分布，以验证测试是否平衡。
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50) #k值的分布不太平衡
# 首先检查未校正 P 值的分布，以验证测试是否平衡。
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50) #k值的分布不太平衡

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)


#绘制knn图
psg_milo <- buildNhoodGraph(psg_milo)
#绘制最后得到的KNN图
plotUMAP(psg_milo) + plotNhoodGraphDA(psg_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")

plotNhoodMA(da_results, alpha = 0.05)

#根据检测到的细胞社区自动聚类
set.seed(29)
da_results <- groupNhoods(psg_milo, da_results, max.lfc.delta = 3, overlap = 2, merge.discord=TRUE)
head(da_results) 
plotNhoodGroups(psg_milo, da_results, layout="UMAP") 

# 使用MATLAB的GenLouvain算法聚类后得到的结果
Nhood_group_data <- read_csv("~/Documents/MATLAB/Nhood_group_updated.csv")
Nhood_group = as.matrix(Nhood_group_data$NhoodGroup_update)
da_results$NhoodGroup <- as.character(Nhood_group)
plotNhoodGroups(psg_milo, da_results, layout="UMAP") 

# 找到每个Nhood的group和cell
nhoods_matrix <- as.matrix(nhoods(psg_milo))

# 假设 dataset 是一个数据框，其中列名为特征名，行名为细胞ID
# 转化 nhoods_matrix 为适合操作的格式
nhoods_matrix <- as.data.frame(nhoods_matrix)

# 转化逻辑矩阵为细胞列表
nhood_cell_indices <- lapply(seq_len(ncol(nhoods_matrix)), function(col_idx) {
  # 转换为逻辑值，假设1为TRUE, 0为FALSE
  logical_column = nhoods_matrix[, col_idx] == 1
  which(logical_column)
})



# 将每个病人的Nhoodgroup信息导出
afib_predict <- psg_milo@colData@listData[["afibPredict"]]

# 1. 计算每个细胞所属的所有 Nhood
nhoods_matrix <- as.matrix(nhoods(psg_milo))
nhoods_matrix <- as.data.frame(nhoods_matrix)
logical_matrix <- nhoods_matrix == 1

# 计算每个细胞所属的所有 Nhood
cell_nhoods <- apply(logical_matrix, 1, function(row) {
  names(which(row))
})

# 2. 根据 nhoodIndex 将 Nhood 转换为 NhoodGroup
nhood_to_group <- setNames(da_results$NhoodGroup, da_results$Nhood)

# 3. 将每个细胞的所有 Nhood 转换为 NhoodGroup，并计算每个细胞的最常见的 NhoodGroup
# 反转 nhoodIndex 映射，从 Nhood 索引映射到 Nhood
nhoodIndex=nhoodIndex(psg_milo)
inverse_nhood_index <- setNames(seq_along(nhoodIndex), unlist(nhoodIndex))
# 使用映射转换 cell_nhoods 中的 Nhood 索引到 Nhood
# converted_nhoods <- inverse_nhood_index[cell_nhoods]
mapped_nhoods <- lapply(cell_nhoods, function(indices) {
  # 过滤掉任何不存在于 inverse_nhood_index 中的索引
  valid_indices <- indices[indices %in% names(inverse_nhood_index)]
  if (length(valid_indices) > 0) {
    # 转换索引到 Nhood 编号
    return(inverse_nhood_index[valid_indices])
  } else {
    return(NA)  # 如果没有有效的索引，返回 NA
  }
})


# 根据每个细胞最邻近的Nhood距离（euclidean距离）
distances_list = list()
for (i in seq_along(cell_nhoods)) {
  cell_coords <- umap_embeddings[i, c("umap_1", "umap_2")]  # 获取当前细胞的UMAP坐标
  
  # 为每个邻域计算距离
  distances <- sapply(cell_nhoods[[i]], function(nhood_index) {
    if (nhood_index %in% rownames(umap_embeddings)) {
      nhood_coords <- umap_embeddings[nhood_index, c("umap_1", "umap_2")]  # 获取邻域中心细胞的UMAP坐标
      # 计算距离
      dist <- sqrt(sum((nhood_coords - cell_coords)^2))
      return(dist)
    } else {
      NA  # 如果没有找到对应的邻域中心细胞，返回NA
    }
  })
  distances_list[[i]]=distances
}

# 初始化一个列表来存储每个细胞的最小距离邻域的NhoodGroup
cell_groups <- vector("character", length = length(cell_nhoods))

# 处理每个细胞
for (i in seq_along(distances_list)) {
  # 获取当前细胞的距离列表并排序
  sorted_indices <- order(distances_list[[i]])
  sorted_nhood_indices <- cell_nhoods[[i]][sorted_indices]
  sorted_nhoods <- inverse_nhood_index[sorted_nhood_indices]
  
  # 初始化当前细胞的NhoodGroup为NA
  cell_groups[i] <- NA
  
  # 检查每个排序后的邻域
  for (nhood in sorted_nhoods) {
    if (any(da_results$Nhood == nhood)) {
      # 检查显著性
      if (da_results$SpatialFDR[da_results$Nhood == nhood] < 0.1) {
        cell_groups[i] <- da_results$NhoodGroup[da_results$Nhood == nhood]
        break  # 找到显著的邻域后停止检查
      }
    }
  }
}

cell_nhood_group <- vector("character", length = length(afib_predict))

# 设置 cell_nhood_group 的值
for (i in seq_along(afib_predict)) {
  if (afib_predict[i] == 1) {
    cell_nhood_group[i] <- cell_groups[i]  # 使用之前确定的组别
  } else {
    cell_nhood_group[i] <- "0"  # afib_predict 不是 1，则组别设置为 0
  }
}

# 计算每个afib=1的病人最近的显著的nhood
# 1. 过滤出显著的 nhood
significant_nhoods <- da_results[da_results$SpatialFDR < 0.05, ]

# 2. 计算距离并找到最近的 nhood
afib_cells <- which(afib_predict == 1)  # 找到afib_predict为1的细胞索引
min_distances <- numeric(length(afib_cells))
min_nhood_indices <- numeric(length(afib_cells))
min_nhood_ids <- numeric(length(afib_cells))
min_nhood_groups  <- numeric(length(afib_cells))

for (i in seq_along(afib_cells)) {
  cell_index <- afib_cells[i]
  # 提取当前细胞的坐标
  cell_coords <- umap_embeddings[cell_index, c("umap_1", "umap_2")]
  
  # 对所有显著nhood计算距离
  distances <- numeric(nrow(significant_nhoods))
  for (j in seq_len(nrow(significant_nhoods))) {
    nhood_index <- as.numeric(significant_nhoods$Nhood[j])
    actual_positions <- nhoodIndex[[nhood_index]]
    nhood_coords <- umap_embeddings[actual_positions, c("umap_1", "umap_2")]
    distances[j] <- sqrt(sum((cell_coords - nhood_coords)^2))
  }
  
  # 找到最小距离及其对应的 nhood 索引和编号
  if (length(distances) > 0) {
    min_index <- which.min(distances) #这是在所有显著的Nhood中的index
    min_distances[i] <- distances[min_index] #这是距离最小值
    min_nhood_indices[i] <- as.numeric(min_index) #这是距离最小值的坐标
    min_nhood_ids[i] <- significant_nhoods$Nhood[min_index]
    min_nhood_groups[i] <- significant_nhoods$NhoodGroup[min_index]
    cat("Processing cell", i)
  } else {
    min_distances[i] <- NA
    min_nhood_indices[i] <- NA
    min_nhood_ids[i] <- NA
    min_nhood_groups[i] <- NA
  }
}

# 3. 结果汇总
results_df <- data.frame(
  CellIndex = afib_cells,
  MinDistance = min_distances,
  MinNhoodIndex = min_nhood_indices,
  MinNhoodID = min_nhood_ids
)

matched_indices <- match(results_df$MinNhoodID, da_results$Nhood)
results_df$MinNhoodGroup <- da_results$NhoodGroup[matched_indices]

# 设置 cell_nhood_group 的值
cell_nhood_group <- vector("character", length = length(afib_predict))

# 使用向量化的 ifelse 函数分配组别
cell_nhood_group <- ifelse(afib_predict == 1, 
                           results_df$MinNhoodGroup[match(seq_along(afib_predict), results_df$CellIndex)],
                           "0")


# 结果赋值给 umap_embeddings 的新列
# 转换为因子
umap_embeddings$NhoodGroup <- cell_nhood_group
# 输出结果



## 计算邻接矩阵
.build_nhood_adjacency <- function(nhoods, overlap=1){
  nh_intersect_mat <- Matrix::crossprod(nhoods)
  nh_intersect_mat[nh_intersect_mat < overlap] <- 0
  
  rownames(nh_intersect_mat) <- colnames(nhoods)
  colnames(nh_intersect_mat) <- colnames(nhoods)
  return(nh_intersect_mat)
}
nhoods_matrix <- nhoods(psg_milo)
nh_intersect_mat <- .build_nhood_adjacency(nhoods_matrix, overlap=1)
nh_intersect_mat<- as.matrix(nh_intersect_mat)


# 将 nhoodIndex 从列表转换为数据框，这里假设 nhoodIndex 是一个命名列表
nhoodIndex=nhoodIndex(psg_milo)
nhoodIndex_df <- data.frame(Nhood = 1:length(nhoodIndex), NhoodIndice = unlist(nhoodIndex), stringsAsFactors = FALSE)
# 合并 da_results 和 nhoodIndex_df
merged_df <- merge(da_results, nhoodIndex_df, by = "Nhood")
# 由于 umap_embeddings 索引与 NhoodIndice 相对应，确保它是一个数据框并且有行名作为 NhoodIndice
umap_embeddings_df <- data.frame(umap_embeddings, NhoodIndice = rownames(umap_embeddings))
# 最终合并所有信息
merged_df <- merge(merged_df, umap_embeddings_df, by = "NhoodIndice")
# 查看最终的合并数据框
head(merged_df)

# 导出数据

