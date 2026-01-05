# 加载包
library(rtracklayer)
library(dplyr)
library(tidyverse)
library(data.table)

# 设置工作目录（根据实际情况修改）
setwd("F:\\ycs\\工作\\work_transcriptome")

# 函数：从GTF文件计算内含子长度
calculate_intron_lengths <- function(gtf_file) {
  # 读取GTF文件
  gtf <- import(gtf_file)
  
  # 转换为data.frame以便处理
  gtf_df <- as.data.frame(gtf)
  
  # 只保留外显子和基因记录
  exons <- gtf_df %>% 
    filter(type == "exon") %>%
    arrange(seqnames, gene_id, transcript_id, start)
  
  # 初始化存储内含子长度的向量
  intron_lengths <- numeric(0)
  
  # 按转录本分组计算内含子
  transcripts <- unique(exons$transcript_id)
  
  for (transcript in transcripts) {
    # 获取该转录本的所有外显子
    transcript_exons <- exons %>% 
      filter(transcript_id == transcript) %>%
      arrange(start)
    
    # 如果该转录本有多个外显子
    if (nrow(transcript_exons) > 1) {
      for (i in 1:(nrow(transcript_exons) - 1)) {
        # 计算内含子长度：下一个外显子的起始 - 当前外显子的结束 - 1
        intron_start <- transcript_exons$end[i] + 1
        intron_end <- transcript_exons$start[i + 1] - 1
        intron_length <- intron_end - intron_start + 1
        
        # 确保内含子长度为正
        if (intron_length > 0) {
          intron_lengths <- c(intron_lengths, intron_length)
        }
      }
    }
  }
  
  return(intron_lengths)
}

# 方法2：更高效的方法（使用data.table）
calculate_intron_lengths_fast <- function(gtf_file) {
  # 读取GTF文件
  gtf <- import(gtf_file)
  gtf_df <- as.data.frame(gtf)
  
  # 转换为data.table
  setDT(gtf_df)
  
  # 筛选外显子并按转录本排序
  exons <- gtf_df[type == "exon"]
  exons <- exons[order(transcript_id, start)]
  
  # 计算内含子长度
  introns <- exons[, {
    if (.N > 1) {
      # 计算当前外显子结束到下一个外显子开始之间的距离
      intron_starts = end[1:(.N-1)] + 1
      intron_ends = start[2:.N] - 1
      intron_lengths = intron_ends - intron_starts + 1
      
      # 只保留正长度的内含子
      valid_introns = intron_lengths > 0
      list(intron_length = intron_lengths[valid_introns])
    }
  }, by = transcript_id]
  
  return(introns$intron_length)
}

# 函数：绘制内含子长度分布
plot_intron_distribution <- function(intron_lengths, 
                                     sample_name = "Sample",
                                     log_transform = TRUE,
                                     binwidth = NULL) {
  
  # 创建数据框
  df <- data.frame(length = intron_lengths)
  
  # 移除极端值（可选）
  # df <- df %>% filter(length <= quantile(length, 0.99))
  
  # 基本统计信息
  cat("内含子长度统计信息:\n")
  cat(paste("样本:", sample_name, "\n"))
  cat(paste("内含子总数:", length(intron_lengths), "\n"))
  cat(paste("平均长度:", round(mean(intron_lengths), 2), "bp\n"))
  cat(paste("中位数:", median(intron_lengths), "bp\n"))
  cat(paste("最小值:", min(intron_lengths), "bp\n"))
  cat(paste("最大值:", max(intron_lengths), "bp\n"))
  cat(paste("标准差:", round(sd(intron_lengths), 2), "bp\n"))
  
  # 创建绘图
  if (log_transform) {
    df$log_length <- log10(df$length + 1)
    
    p <- ggplot(df, aes(x = log_length)) +
      geom_density(fill = "steelblue", alpha = 0.7) +
      labs(
        title = paste("内含子长度分布密度图 -", sample_name),
        subtitle = "X轴为log10(长度+1)",
        x = "log10(内含子长度 + 1)",
        y = "密度"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.title = element_text(size = 12)
      )
  } else {
    p <- ggplot(df, aes(x = length)) +
      geom_density(fill = "coral", alpha = 0.7) +
      labs(
        title = paste("内含子长度分布密度图 -", sample_name),
        x = "内含子长度 (bp)",
        y = "密度"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12)
      ) +
      scale_x_continuous(labels = scales::comma)
  }
  
  # 添加直方图叠加（可选）
  p <- p + geom_histogram(aes(y = ..density..), 
                          bins = ifelse(is.null(binwidth), 50, NULL),
                          binwidth = binwidth,
                          alpha = 0.3, 
                          fill = "gray50")
  
  return(p)
}

# 主程序
main <- function() {
  # 设置GTF文件路径
  gtf_file <- "Gmax_ZH13_v2.0.gene.gtf"  # 替换为您的GTF文件路径
  
  # 检查文件是否存在
  if (!file.exists(gtf_file)) {
    stop("GTF文件不存在，请检查文件路径")
  }
  
  cat("正在读取GTF文件并计算内含子长度...\n")
  
  # 使用方法1（适用于小型GTF文件）
  # intron_lengths <- calculate_intron_lengths(gtf_file)
  
  # 使用方法2（更高效，适用于大型GTF文件）
  intron_lengths <- calculate_intron_lengths_fast(gtf_file)
  
  cat(paste("成功计算", length(intron_lengths), "个内含子长度\n"))
  
  # 绘制分布图
  cat("正在绘制分布图...\n")
  
  # 绘制对数转换的分布图
  p1 <- plot_intron_distribution(intron_lengths, 
                                 sample_name = basename(gtf_file),
                                 log_transform = TRUE)
  
  # 绘制原始尺度分布图（使用分位数限制显示范围）
  df <- data.frame(length = intron_lengths)
  q99 <- quantile(intron_lengths, 0.99)
  df_filtered <- df %>% filter(length <= q99)
  
  p2 <- ggplot(df_filtered, aes(x = length)) +
    geom_density(fill = "lightgreen", alpha = 0.7) +
    labs(
      title = paste("内含子长度分布（99%分位数内） -", basename(gtf_file)),
      x = "内含子长度 (bp)",
      y = "密度"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    ) +
    scale_x_continuous(labels = scales::comma)
  
  # 保存图形
  ggsave("intron_length_distribution_log.pdf", p1, width = 10, height = 6)
  ggsave("intron_length_distribution_raw.pdf", p2, width = 10, height = 6)
  
  cat("图形已保存为PDF文件\n")
  
  # 显示图形
  print(p1)
  print(p2)
  
  # 保存内含子长度数据
  write.csv(data.frame(intron_length = intron_lengths),
            "intron_lengths.csv", 
            row.names = FALSE)
  
  cat("内含子长度数据已保存到intron_lengths.csv\n")
}

# 运行主程序
main()

# 可选：创建更详细的统计摘要
create_summary_statistics <- function(intron_lengths) {
  summary_df <- data.frame(
    Statistic = c("总数", "平均值", "中位数", "最小值", "最大值", 
                  "标准差", "第一分位数", "第三分位数"),
    Value = c(
      length(intron_lengths),
      round(mean(intron_lengths), 2),
      median(intron_lengths),
      min(intron_lengths),
      max(intron_lengths),
      round(sd(intron_lengths), 2),
      quantile(intron_lengths, 0.25),
      quantile(intron_lengths, 0.75)
    )
  )
  
  return(summary_df)
}

# 如果只需要快速查看分布，可以使用这个简化版本
quick_plot <- function(gtf_file) {
  # 快速计算内含子长度（只处理前1000个转录本以加快速度）
  gtf <- import(gtf_file)
  gtf_df <- as.data.frame(gtf)
  exons <- gtf_df %>% filter(type == "exon")
  
  # 获取前1000个转录本
  transcripts <- unique(exons$transcript_id)[1:min(1000, length(unique(exons$transcript_id)))]
  exons_sample <- exons %>% filter(transcript_id %in% transcripts)
  
  intron_lengths <- numeric(0)
  
  for (transcript in transcripts) {
    transcript_exons <- exons_sample %>% 
      filter(transcript_id == transcript) %>%
      arrange(start)
    
    if (nrow(transcript_exons) > 1) {
      for (i in 1:(nrow(transcript_exons) - 1)) {
        intron_length <- transcript_exons$start[i + 1] - transcript_exons$end[i] - 1
        if (intron_length > 0) {
          intron_lengths <- c(intron_lengths, intron_length)
        }
      }
    }
  }
  
  # 快速绘图
  hist(log10(intron_lengths + 1), 
       main = "内含子长度分布 (log10)",
       xlab = "log10(长度+1)",
       col = "lightblue",
       breaks = 50)
  
  return(intron_lengths)
}