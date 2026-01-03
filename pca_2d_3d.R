library("tidyverse")
pca_data <- read.table("F:\\ycs\\select_cluster\\sv_pca\\pca_results.eigenvec", header = FALSE)
head(pca_data)
p1 <- ggplot(pca_data, aes(x = pca_data$V3, y = pca_data$V4)) +
  geom_point() +
  theme_classic()
ggsave("F:\\ycs\\select_cluster\\sv_pca\\pca_results_2D.eigenvec.pdf", p1, width = 8, height = 8)
library(scatterplot3d)
pdf("F:\\ycs\\select_cluster\\sv_pca\\pca_results_3D.eigenvec.pdf", width = 8, height = 8)
p2 <- scatterplot3d(x = pca_data$V3,
                     y = pca_data$V4,
                     z = pca_data$V5,
                     main = "3D PCA Scatter Plot with Labels",
                     xlab = "PC1 (V3)",
                     ylab = "PC2 (V4)",
                     zlab = "PC3 (V5)",
                     pch = 16,
                     color = "black",
                     angle = 55,
                     box = TRUE)  # 直接在这里设置为TRUE

# 添加文本标签
text(p2$xyz.convert(pca_data$V3, pca_data$V4, pca_data$V5),
     labels = pca_data$V1,
     cex = 1,  # 字体大小
     pos = 4,    # 位置：4=右侧
     col = adjustcolor("black", alpha.f = 0.4))
dev.off()



