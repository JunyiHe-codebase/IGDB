pca_data <- read.table("F:\\ycs\\select_cluster\\sv_pca\\pca_results_236.eigenvec", header = FALSE)
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
head(result)
write.table(result, "F:\\ycs\\select_cluster\\sv_pca\\pca_results_236.eigenvec.cluster.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
library("tidyverse")
library("scatterplot3d")
p1 <- ggplot(result, aes(x = pca_1, y = pca_2, color = as.factor(Cluster))) +
  geom_point() +
  theme_classic()
ggsave("F:\\ycs\\select_cluster\\sv_pca\\pca_results_236.eigenvec.2d.pdf", p1, width = 8, height = 8)
pdf("F:\\ycs\\select_cluster\\sv_pca\\pca_results_236.eigenvec.3d.pdf", width = 8, height = 8)
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
dev.off()
