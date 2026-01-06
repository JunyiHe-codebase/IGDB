library(tidyverse)
library(dplyr)
library(tidyr)
library(gridExtra)

# 创建数据框
data <- read.table("F:\\ycs\\select_cluster\\sv_pca\\30_sample_maf005_geno09_filtered.txt", header = T, sep = "\t")

# 确保Sample的顺序保持不变
data$Sample <- factor(data$Sample, levels = data$Sample)

# 1. 绘制SVs的柱状图
plot1 <- ggplot(data, aes(x = Sample, y = SVs)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  ) +
  labs(title = "SVs per Sample", x = "Sample", y = "SVs") +
  scale_y_continuous(labels = scales::comma)

# 2. 绘制num_SVs_captured的柱状图
plot2 <- ggplot(data, aes(x = Sample, y = num_SVs_captured)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  ) +
  labs(title = "Cumulative SVs Captured", x = "Sample", y = "Cumulative SVs") +
  scale_y_continuous(labels = scales::comma)

# 3. 绘制SVs_captured的柱状图
plot3 <- ggplot(data, aes(x = Sample, y = SVs_captured)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  ) +
  labs(title = "Proportion of SVs Captured", x = "Sample", y = "Proportion") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))

# 将三个图形组合在一起
grid.arrange(plot1, plot2, plot3, ncol = 3)

# 如果想要更大的图形并保存为文件，可以使用以下代码：
# 设置图形大小
options(repr.plot.width = 18, repr.plot.height = 6)

# 保存图形
ggsave("F:\\ycs\\select_cluster\\sv_pca\\30_sample_maf005_geno09_filtered.pdf", 
       arrangeGrob(plot1, plot2, plot3, ncol = 3),
       width = 18, height = 6)