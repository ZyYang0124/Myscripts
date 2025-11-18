# ================================
# R 脚本：基于 Phyluce 统计结果生成 UCE 数量分布小提琴图
# 作者：Qwen
# ================================

# 加载必要的包
library(ggplot2)
library(dplyr)
library(readr)

# 设置输入文件路径
input_file <- "C:/Users/l/Desktop/313spp.UCE.statistics.csv"

# 自动获取输出目录（与输入文件同位置）
output_dir <- dirname(input_file)

# 定义列名
col_names <- c("filename", "n", "total_length", "mean", "std_dev", 
               "shortest", "longest", "median", "N50")

# 读取数据
df <- read_delim(input_file, delim = ",", col_names = col_names, show_col_types = FALSE)

# 提取物种名
df$species <- tools::file_path_sans_ext(basename(df$filename))

# 提取属名（假设物种名格式为 'Genus-species'）
df$genus <- sapply(strsplit(df$species, "-"), `[`, 1)  # 按 '-' 分割字符串，取第一部分作为属名

# 计算统计值
mean_n <- round(mean(df$n), 1)
median_n <- median(df$n)
min_n <- min(df$n)
max_n <- max(df$n)
sd_n <- round(sd(df$n), 1)

# 绘制按属分组的小提琴图
p <- ggplot(df, aes(x = genus, y = n)) +  # 使用 genus 作为 x 轴
  geom_violin(fill = "lightgray", color = "black", alpha = 0.8,
              draw_quantiles = c(0.25, 0.5, 0.75)) +  # 添加四分位线
  geom_jitter(width = 0.1, alpha = 0.6, color = "darkred", size = 1.2) +  # 数据点
  geom_hline(yintercept = mean_n, color = "blue", linetype = "dashed", size = 1) +
  annotate("text", x = length(unique(df$genus))/2, y = mean_n + (max_n - min_n)*0.02,
           label = paste("Mean =", mean_n), color = "blue", size = 4) +
  labs(
    title = "Distribution of UCE Loci Count by Genus",
    subtitle = paste("Range:", min_n, "–", max_n, "| Median =", median_n, "| SD =", sd_n),
    y = "Number of UCE Loci",
    x = "Genus"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转x轴标签以避免重叠
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  )

# 显示图形
print(p)

# 保存图表
png_output <- file.path(output_dir, "uce_count_distribution_by_genus_violin.png")
pdf_output <- file.path(output_dir, "uce_count_distribution_by_genus_violin.pdf")

ggsave(png_output, p, width = 40, height = 6, dpi = 300, bg = "white")  # 调整宽度适应更多的属名
ggsave(pdf_output, p, width = 40, height = 6, bg = "white")

cat("✅ 小提琴图 (PNG) 已保存到：", png_output, "\n")
cat("✅ 小提琴图 (PDF) 已保存到：", pdf_output, "\n")

# 输出统计摘要
cat("\n📊 UCE 数量分布统计：\n")
cat("  总物种数: ", nrow(df), "\n")
cat("  范围: ", min_n, " – ", max_n, "\n")
cat("  中位数: ", median_n, "\n")
cat("  平均值: ", mean_n, " ± ", sd_n, "\n")