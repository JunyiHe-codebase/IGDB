library(tidyverse)
library(dplyr)
library(tidyr)

df <- read.csv("F:\\ycs\\select_cluster\\2025.csv", header = T)
cluster <- read.table("F:\\ycs\\select_cluster\\cluster.txt", header = T, sep = "\t")
cluster_sample <- cluster$IID
extract_data <- df %>% filter(ID %in% cluster_sample)
head(extract_data)

#处理
deal_data <- extract_data[c("ID", "Seed_coat_color")]
merge_data <- merge(deal_data, cluster, 
                                by.x = names(deal_data)[1], 
                                by.y = names(cluster)[1])
merge_data_clean <- merge_data %>% 
  filter(!is.na(Seed_coat_color) & !is.na(Cluster) & !is.na(is_special))
head(merge_data_clean)

merge_data_clean <- merge_data_clean %>%
  mutate(
    # 将Cluster转换为因子
    Cluster = factor(Cluster),
    # 确保Seed_coat_color是数值型
    Seed_coat_color = factor(Seed_coat_color)
  )


p1 <- ggplot(merge_data_clean, aes(x = Cluster, fill = Seed_coat_color)) +
  geom_bar(aes(group = Seed_coat_color), 
           position = position_dodge2(preserve = "single", padding = 0.1)) +
  # 添加计数文本
  geom_text(
    aes(group = Seed_coat_color, label = after_stat(count)),
    stat = "count",
    position = position_dodge2(width = 0.9, preserve = "single"),
    vjust = -0.5,  # 将文本放在柱子上方
    size = 3
  ) +
  # 调整y轴范围，为文本留出空间
  ylim(0, max(table(merge_data_clean$Cluster, merge_data_clean$Seed_coat_color)) * 1.1)
ggsave("F:\\ycs\\select_cluster\\Seed_coat_color.pdf", p1, width = 6, height = 10)
