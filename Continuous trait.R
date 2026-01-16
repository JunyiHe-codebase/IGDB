library(dplyr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(coin)
library(patchwork)

df <- read.csv("F:\\ycs\\select_cluster\\2025.csv", header = T)
cluster <- read.table("F:\\ycs\\select_cluster\\cluster.txt", header = T, sep = "\t")
cluster_sample <- cluster$IID
extract_data <- df %>% filter(ID %in% cluster_sample)
head(extract_data)

#处理指定数据列
deal_data <- extract_data[c("ID", "grain_weight")]
merge_data <- merge(deal_data, cluster, 
                                by.x = names(deal_data)[1], 
                                by.y = names(cluster)[1])
merge_data_clean <- merge_data %>% 
  filter(!is.na(grain_weight) & !is.na(Cluster) & !is.na(is_special))
set.seed(123)  # 设置随机种子以保证jitter的可重复性
# 创建主要图形（使用非参数检验方法）


if(length(unique(merge_data_clean$Cluster)) == 2) {
  # 执行Mann-Whitney U检验
  mw_test <- wilcox.test(grain_weight ~ Cluster, data = merge_data_clean, 
                         exact = FALSE)# 对于大样本使用近似计算
  # 计算效应量 (r)
  effect_size <- wilcox_effsize(merge_data_clean, grain_weight ~ Cluster)
  effect_size_label <- sprintf("r = %.3f", effect_size$effsize)

  # 提取p值用于图形标注
  p_value_mw <- mw_test$p_value
  p_label_mw <- ifelse(p_value_mw < 0.001, "p < 0.001", 
                    ifelse(p_value_mw < 0.01, sprintf("p = %.3f", p_value_mw),
                           sprintf("p = %.3f", p_value_mw)))
}

# 计算y轴范围用于显著性标记
y_max <- max(merge_data_clean$grain_weight, na.rm = TRUE)
y_min <- min(merge_data_clean$grain_weight, na.rm = TRUE)
y_range <- y_max - y_min
y_pos_line <- y_max + 0.12 * y_range
y_pos_text <- y_max + 0.2 * y_range

# 创建基础箱型图
p1 <- ggplot(merge_data_clean, aes(x = factor(Cluster), y = grain_weight, 
                          fill = factor(Cluster))) +
  # 添加箱型图
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA, 
               color = "black", linewidth = 0.5) +
  
  # 添加散点（所有样本点）
  geom_jitter(aes(color = is_special == "Target", shape = is_special == "Target"), 
              width = 0.2, height = 0, size = 2.5, alpha = 0.8) +
  
  # 设置颜色：Target样本点为红色，其他为灰色
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray40"),
                     name = "Sample Type",
                     labels = c("TRUE" = "Target", "FALSE" = "Other"),
                     guide = guide_legend(override.aes = list(size = 3))) +
  
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
                     name = "Sample Type",
                     labels = c("TRUE" = "Target", "FALSE" = "Other")) +
  
  scale_fill_manual(values = c("1" = "#66C2A5", "2" = "#FC8D62"),
                    name = "Cluster") +
  
  # 添加比较线和显著性标记
  geom_segment(x = 1, xend = 2, y = y_pos_line, yend = y_pos_line, 
               color = "black", linewidth = 0.8) +
  geom_segment(x = 1, xend = 1, y = y_pos_line - 0.01*y_range, yend = y_pos_line, 
               color = "black", linewidth = 0.8) +
  geom_segment(x = 2, xend = 2, y = y_pos_line - 0.01*y_range, yend = y_pos_line, 
               color = "black", linewidth = 0.8) +
  
  annotate("text", x = 1.5, y = y_pos_text, 
           label = paste0("Mann-Whitney U Test\n", p_label_mw),
           size = 4, fontface = "bold", color = ifelse(p_value_mw < 0.05, "red", "black")) +
  
  annotate("text", x = 1.5, y = y_pos_text - 0.08 * y_range, 
           label = effect_size_label,
           size = 3.5, color = "darkgreen") +
  
  # 设置主题和标签
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.2),
    axis.text = element_text(color = "black")
  ) +
  
  labs(
    title = "Non-parametric Test (Mann-Whitney U)",
    x = "Cluster",
    y = "grain_weight"
  ) +
  
  # 确保所有图形y轴范围一致
  ylim(y_min, y_pos_text + 0.05 * y_range)

# 创建第二个图形(p2) - 参数检验版本
# 执行t检验
if(length(unique(merge_data_clean$Cluster)) == 2) {
  # 方差齐性检验
  var_test <- var.test(grain_weight ~ Cluster, data = merge_data_clean)
  var_equal <- var_test$p.value >= 0.05
  
  # 执行t检验
  t_test <- t.test(grain_weight ~ Cluster, 
                   data = merge_data_clean, 
                   var.equal = var_equal)
  p_value_t <- t_test$p.value
  
  # 创建p2图形（t检验版本）
  p2 <- ggplot(merge_data_clean, 
               aes(x = factor(Cluster), y = grain_weight, 
                   fill = factor(Cluster))) +
    # 添加箱型图
    geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA, 
                 color = "black", linewidth = 0.5) +
    
    # 添加所有样本点
    geom_jitter(aes(color = is_special == "Target", shape = is_special == "Target"), 
                width = 0.2, height = 0, size = 2.5, alpha = 0.8) +
    
    # 设置颜色和形状
    scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "gray50"),
                       name = "Sample Type",
                       labels = c("TRUE" = "Target", "FALSE" = "Other")) +
    
    scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
                       name = "Sample Type",
                       labels = c("TRUE" = "Target", "FALSE" = "Other")) +
    
    scale_fill_manual(values = c("1" = "#66C2A5", "2" = "#FC8D62"),
                      name = "Cluster") +
    
    # 添加比较线和显著性标记
    geom_segment(x = 1, xend = 2, y = y_pos_line, yend = y_pos_line, 
                 color = "black", linewidth = 0.8) +
    geom_segment(x = 1, xend = 1, y = y_pos_line - 0.01*y_range, yend = y_pos_line, 
                 color = "black", linewidth = 0.8) +
    geom_segment(x = 2, xend = 2, y = y_pos_line - 0.01*y_range, yend = y_pos_line, 
                 color = "black", linewidth = 0.8) +
    
    # 添加t检验结果
    annotate("text", x = 1.5, y = y_pos_text, 
             label = paste0(ifelse(var_equal, "Student's t-test\n", "Welch's t-test\n"), 
                            ifelse(p_value_t < 0.001, "p < 0.001", 
                                   sprintf("p = %.3f", p_value_t))),
             size = 4, fontface = "bold", color = ifelse(p_value_t < 0.05, "blue", "black")) +
    
    annotate("text", x = 1.5, y = y_pos_text - 0.08 * y_range, 
             label = ifelse(var_equal, "Equal variances", "Unequal variances"),
             size = 3.5, color = "darkorange") +
    
    # 设置主题和标签
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.2),
      axis.text = element_text(color = "black"),
      axis.title.y = element_blank()  # 移除y轴标题，因为两个图共享
    ) +
    
    labs(
      title = "Parametric Test (t-test)",
      x = "Cluster"
    ) +
    
    # 确保所有图形y轴范围一致
    ylim(y_min, y_pos_text + 0.05 * y_range)
}

# 组合图形 - 并排排列
combined_plot <- p1 + p2 + 
  plot_layout(ncol = 2) +  # 两列布局
  plot_annotation(
    title = "Comparison of grain_weight Between Clusters",
    subtitle = "Target samples are highlighted as red triangles",
    tag_levels = 'A',  # 添加A/B标签
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 12))
  )

# 添加图例（只显示一次）
legend_plot <- cowplot::get_legend(
  p1 + 
    theme(legend.position = "bottom",
          legend.box = "horizontal") +
    guides(color = guide_legend(title.position = "top", nrow = 1),
           shape = guide_legend(title.position = "top", nrow = 1),
           fill = guide_legend(title.position = "top", nrow = 1))
)

# 最终组合图形
final_plot <- combined_plot / legend_plot + 
  plot_layout(heights = c(10, 1))

# 显示图形
print(final_plot)

# 保存图形
ggsave("F:\\ycs\\select_cluster\\grain_weight.pdf", 
       plot = final_plot,
       width = 16, 
       height = 9,
       bg = "white")

# 输出统计结果
cat("=== 统计检验结果汇总 ===\n")
cat(sprintf("Mann-Whitney U 检验: p = %.4f\n", p_value_mw))
cat(sprintf("效应量 (r): %.3f\n", effect_size$effsize))
cat(sprintf("t检验: p = %.4f\n", p_value_t))
cat(sprintf("方差齐性: %s (p = %.4f)\n", 
            ifelse(var_equal, "是", "否"), var_test$p.value))
cat(sprintf("推荐检验方法: %s\n", 
            ifelse(shapiro.test(merge_data_clean$grain_weight)$p.value < 0.05, 
                   "非参数检验", "参数检验")))
