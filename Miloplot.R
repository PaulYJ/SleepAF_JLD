# 导出所有的图片：
install.packages("ggbeeswarm")
library("ggbeeswarm")
install.packages("dplyr")
library(dplyr)

# milo DA on UMAP visualizaion
p <- plotUMAP(psg_milo) + 
  plotNhoodGraphDA(psg_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")

# milo DA bee swarm
p <- plotDAbeeswarm(da_results, "NhoodGroup",alpha=0.05) #绘制不同cluster的log fold change

# milo Nhoodsize
p <- plotNhoodSizeHist(psg_milo)

# milo Nhood cluster
plot_groups <- unique(da_results$NhoodGroup)
groups_res <- da_results[da_results$NhoodGroup %in% plot_groups,]
colData(psg_milo)[unlist(nhoodIndex(psg_milo)[groups_res$Nhood]),"NhoodGroup"] <- groups_res$NhoodGroup
# 定义颜色映射，假设有6个正常组别加上NaN="0"
colors <- c("NaN" = "gray", "1" = "#E64C4C", "2" = "#E6B422", "3" = "#4DB24C", "4" = "#1A99CC", "5" = "#333399", "6" = "#B24C98")
library(ggraph)
library(igraph)
plot <- plotNhoodGraph(psg_milo, colour_by = "NhoodGroup")
p <-plot + scale_fill_manual(values = colors)

# sleep parameters beeswarm
# 变量设置和描述性标题
vars_settings <- list(
  slpprdp = c('Total Sleep Time (min)', 10, 540),
  waso = c('WASO (min)', 10, 540),
  slplatp = c('Latency (min)', 5, 200),
  slpeffp = c('Efficiency (%)', 3, 100),
  timerem = c('REM Duration (min)', 3, 160),
  timest1 = c('N1 Duration (min)', 2, 140),
  timest2 = c('N2 Duration (min)', 5, 400),
  timest34 = c('N3/4 Duration (min)', 4, 250),
  timeremp = c('REM Proportion (%)', 2, 100),
  timest1p = c('N1 Proportion (%)', 2, 100),
  timest2p = c('N2 Proportion (%)', 2, 100),
  times34p = c('N3/4 Proportion (%)', 2, 100),
  hslptawp = c('Wake per hour (/h)', 0.5, 22),
  ai_all = c('Arousal Index (/h)', 2, 80),
  ai_nrem = c('NREM Arousal Index (/h)', 2, 100),
  ai_rem = c('REM Arousal Index (/h)', 2, 80),
  remlaiip = c('REM Latency (min)', 4, 360)
)

results <- list() 

# 遍历每个变量来计算平均值并绘图
for (feature in names(vars_settings)) {
  settings <- vars_settings[[feature]]
  descriptive_title <- settings[1]  # 描述性标题
  
  feature_row_index <- which(rownames(dataset) == feature)
  nhood_feature_avgs <- data.frame(Nhood_indice = colnames(nhoods_matrix))
  
  # 为当前特征计算每个邻域的平均值
  feature_avgs <- sapply(nhood_cell_indices, function(cell_indices) {
    if (length(cell_indices) > 0) {
      mean(as.numeric(dataset[feature_row_index, cell_indices]), na.rm = TRUE)
    } else {
      NA  # 如果邻域中没有细胞，则返回 NA
    }
  })
  
  # 将计算结果添加到结果数据框
  nhood_feature_avgs[[feature]] <- feature_avgs
  
  
  # 假设 da_results 包含每个邻域的分组信息
  nhoodIndex = nhoodIndex(psg_milo)
  plot_data <- data.frame(
    Nhood = seq_along(nhoodIndex),
    Nhood_indice = unlist(nhoodIndex)
  )
  
  # 合并 plot_data 与 da_results，通过 Nhood 连接
  plot_data <- merge(plot_data, nhood_feature_avgs, by = "Nhood_indice")
  plot_data <- merge(plot_data, da_results, by = "Nhood")
  
  # 过滤掉 NaN 组
  plot_data <- plot_data[plot_data$NhoodGroup != "NaN", ]
  
  # 统计每个 NhoodGroup 的 mean, SD 和 N
  group_stats <- plot_data %>%
    group_by(NhoodGroup) %>%
    summarise(
      group_mean = mean(!!sym(feature), na.rm = TRUE),
      group_sd = sd(!!sym(feature), na.rm = TRUE),
      group_n = sum(!is.na(!!sym(feature)))
    )
  
  results[[feature]] <- group_stats
  
  p <- ggplot(plot_data, aes(x = NhoodGroup, y = !!sym(feature), color = NhoodGroup)) +
    geom_quasirandom() +
    labs(x = "Cluster", y = descriptive_title) +
    scale_color_manual(values = colors) +
    theme_classic() +  # 使用经典主题
    theme(
      legend.position = "none",  # 关闭图例
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 0),  # X轴刻度标签不旋转
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(color = "black")
    ) +
    coord_flip()
  
  # 保存图像到文件
  ggsave(filename = paste0(output_file, feature, "_beeswarm_plot.png"), plot = p, width = 10, height = 5, dpi = 300)
  ggsave(filename = paste0(output_file, feature, "_beeswarm_plot.pdf"),plot = p, device = "pdf", width = 6, height = 2.5)
}



# cluster allocation 散点图
colors <- c(
  "0" = "#D3D3D3",  # 更浅的灰色
  "1" = "#E64C4C",
  "2" = "#E6B422",
  "3" = "#4DB24C",
  "4" = "#1A99CC",
  "5" = "#333399",
  "6" = "#B24C98"
)

names(colors) <- as.character(unique_groups)
# 将NhoodGroup转为字符以便于排序
umap_embeddings$NhoodGroup <- as.character(umap_embeddings$NhoodGroup)
umap_embeddings <- umap_embeddings[order(umap_embeddings$NhoodGroup != "0"), ]

library(ggplot2)
p <-ggplot(umap_embeddings, aes(x = umap_1, y = umap_2, color = as.factor(NhoodGroup))) +
  geom_point(alpha = 1, size = 0.5) +
  scale_color_manual(values = colors) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Nhood Group"  # 设置图例标题
  ) + 
  theme_void() +
  theme(
    legend.position = "right",  # 根据需要调整图例位置
    legend.title = element_text(size = 12),  # 调整图例标题的文字大小
    panel.background = element_rect(fill = NA, colour = NA)  # 设置透明背景
  )


# 绘制clustersize柱状图
# 计算每个NhoodGroup的大小，排除NA
plot_data <- umap_embeddings[umap_embeddings$NhoodGroup != "0", ]
cluster_sizes <- plot_data %>%
  group_by(NhoodGroup) %>%
  summarise(ClusterSize = n(), .groups = 'drop')

# 颜色定义
colors <- c(
  "1" = "#E64C4C",
  "2" = "#E6B422",
  "3" = "#4DB24C",
  "4" = "#1A99CC",
  "5" = "#333399",
  "6" = "#B24C98"
)

# 确保NhoodGroup为因子，并且顺序正确
cluster_sizes$NhoodGroup <- factor(cluster_sizes$NhoodGroup, levels = names(colors))

# 绘制柱状图
p <- ggplot(cluster_sizes, aes(x = NhoodGroup, y = ClusterSize, fill = NhoodGroup)) +
  geom_bar(stat = "identity", width = 0.7) +  # 使用固定宽度的柱体
  geom_text(aes(label = ClusterSize), vjust = -0.3, color = colors, size = 3.5) +  # 在柱子上显示数值
  scale_fill_manual(values = colors) +        # 应用自定义颜色
  labs(x = "Nhood Group", y = "Cluster Size") +  # 标签
  scale_y_continuous(limits = c(0, NA)) +     # 设置Y轴开始于0
  theme_classic() +                           # 使用经典主题
  theme(
    axis.text.x = element_text(angle = 0),    # X轴标签不旋转
    legend.position = "none",                 # 不显示图例
    plot.title = element_text(hjust = 0.5)    # 标题居中
  )

# 显示图形
print(p)

