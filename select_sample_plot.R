library("tidyverse")
library("scatterplot3d")
library("ggrepel")
sample_file <- read.table("F:\\ycs\\select_cluster\\sv_pca\\30_sample_maf005_geno09_filtered.txt", header = TRUE)
pca_data <- read.table("F:\\ycs\\select_cluster\\sv_pca\\pca_results_236.eigenvec", header = FALSE)

sample_list <- sample_file$Sample

pca_matrix <- as.matrix(pca_data[,3:4])
set.seed(123)
kmeans_result <- kmeans(pca_matrix, centers = 2)
result <- data.frame(
  IID = pca_data$V2,
  pca_1 = pca_data$V3,
  pca_2 = pca_data$V4,
  pca_3 = pca_data$V5,
  Cluster = kmeans_result$cluster
)
result <- result %>%
  mutate(is_special = ifelse(IID %in% sample_list, "Target", "Other"))

head(result)
write.table(result, "F:\\ycs\\select_cluster\\sv_pca\\pca_results_236.eigenvec.cluster.selected.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

p1 <- ggplot(result, aes(x = pca_1, y = pca_2, color = as.factor(Cluster))) +
  geom_point(aes(shape = is_special, size = is_special), alpha = 0.7) +
  geom_point(data = subset(result, is_special == "Target"), 
             aes(x = pca_1, y = pca_2),
             shape = 17,  # 三角形
             size = 3,    # 更大的点
             color = "red", # 红色
             stroke = 1.5) + # 边框粗细
  # 设置形状和大小比例
  geom_text_repel(data = subset(result, is_special == "Target"),
                  aes(label = IID),
                  color = "black",
                  size = 3,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  max.overlaps = 20) +
  scale_shape_manual(values = c(Target = 17, Other = 16)) +
  scale_size_manual(values = c(Target = 3, Other = 2)) +
  labs(title = "PCA Plot with Special Samples Highlighted",
       x = "PC1",
       y = "PC2",
       color = "Cluster",
       shape = "Sample Type") +
  theme(legend.position = "right")
ggsave("F:\\ycs\\select_cluster\\sv_pca\\pca_results_236.eigenvec.2d.selected.pdf", p1, width = 8, height = 8)
pdf("F:\\ycs\\select_cluster\\sv_pca\\pca_results_236.eigenvec.3d.selected.pdf", width = 8, height = 8)
p2 <- scatterplot3d(x = result$pca_1,
                    y = result$pca_2,
                    z = result$pca_3,
                    main = "3D PCA Scatter Plot with Labels",
                    xlab = "PC1 (V3)",
                    ylab = "PC2 (V4)",
                    zlab = "PC3 (V5)",
                    pch = 16,
                    color = as.factor(result$Cluster),
                    angle = 55,
                    box = TRUE)
selected_idx <- which(result$IID %in% sample_list)
# 为选定样本添加标签
p2$points3d(x = result$pca_1[selected_idx],
            y = result$pca_2[selected_idx],
            z = result$pca_3[selected_idx],
            col = "red",  # 改变颜色以突出显示
            pch = 17,     # 三角形标记
            cex = 1.5)    # 增大点的大小

# 添加文本标签（稍微偏移以避免重叠）
text(p2$xyz.convert(result$pca_1[selected_idx],
                    result$pca_2[selected_idx],
                    result$pca_3[selected_idx]),
     labels = result$IID[selected_idx],
     cex = 0.7,           # 字体大小
     col = "red",         # 标签颜色
     pos = 4)             # 文本位置（4=右侧）

dev.off()