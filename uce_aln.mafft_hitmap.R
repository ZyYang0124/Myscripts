# ================================
# ✅ 修复版 R 脚本：从 MAFFT 比对文件夹构建存在/缺失矩阵 + 热图
# 无需 phyluce，纯 R 实现
# ================================

# 第一步：加载包（首次运行需安装）
# install.packages(c("ggplot2", "dplyr", "tidyr", "ape"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(ape)

# 第二步：设置比对文件夹路径（修改为你的实际路径）
alignment_dir <- "E:/Myworks/Paper_14_Chrysillini_big/mafft.alignments"  # ✅ 路径正确

# 检查目录是否存在
if (!dir.exists(alignment_dir)) {
  stop("❌ 目录不存在！请检查路径是否正确。")
}

# 获取所有 fasta/fa 文件
files <- list.files(alignment_dir, pattern = "\\.fasta$|\\.fa$", full.names = TRUE)
if (length(files) == 0) {
  stop("❌ 在指定目录中未找到 .fasta 或 .fa 文件，请确认文件扩展名是否匹配。")
}

locus_names <- tools::file_path_sans_ext(basename(files))

# 第三步：读取每个比对文件，提取物种名
cat("🔍 正在扫描比对文件...\n")
presence_list <- lapply(seq_along(files), function(i) {
  file <- files[i]
  locus <- locus_names[i]
  
  # 读取 fasta
  seqs <- try(read.dna(file, format = "fasta", as.character = TRUE), silent = TRUE)
  
  if (inherits(seqs, "try-error") || is.null(seqs) || !is.matrix(seqs)) {
    warning(paste("⚠️ 无法读取或格式错误，跳过文件:", file))
    return(data.frame(locus = locus, taxon = character()))
  }
  
  taxa <- rownames(seqs)
  data.frame(locus = locus, taxon = taxa, stringsAsFactors = FALSE)
})

# 过滤掉空结果，再合并
presence_list <- presence_list[vapply(presence_list, nrow, integer(1)) > 0]
if (length(presence_list) == 0) {
  stop("❌ 所有文件都未能成功读取，请检查文件格式是否为标准 FASTA。")
}

presence_df <- do.call(rbind, presence_list)

# 第四步：构建存在/缺失矩阵（1=存在，0=缺失）
cat("📊 正在构建存在/缺失矩阵...\n")
mat_wide <- presence_df %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = locus, values_from = present, values_fill = 0)

# 第五步：排序（按完整性从高到低）
# 排除第一列 'taxon'
mat_data <- mat_wide[, -1]  # 只保留数值列

# 计算每行（物种）的完整性
row_completeness <- rowSums(mat_data)
# 计算每列（位点）的覆盖度
col_completeness <- colSums(mat_data)

# 添加完整性用于排序
mat_wide$row_completeness <- row_completeness

# 排序列名（按列覆盖度降序）
ordered_loci <- names(mat_data)[order(col_completeness, decreasing = TRUE)]

# 排序行（物种按完整性降序），排序列（位点按覆盖度降序）
mat_sorted <- mat_wide %>%
  arrange(desc(row_completeness)) %>%
  select(taxon, all_of(ordered_loci))  # ✅ 正确语法

# 移除临时列
mat_sorted$row_completeness <- NULL

# 第六步：转为长格式绘图
mat_long <- mat_sorted %>%
  pivot_longer(-taxon, names_to = "locus", values_to = "present")

# 第七步：绘制热图
cat("🎨 正在绘制热图...\n")
p <- ggplot(mat_long, aes(x = locus, y = taxon, fill = present)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_gradient(
    low = "white", high = "red",
    limits = c(0, 1),
    na.value = "gray80",
    breaks = c(0, 1), labels = c("Absent", "Present"),
    guide = guide_legend(title = "Status")
  ) +
  labs(
    title = "UCE Data Matrix Completeness",
    subtitle = paste(nrow(mat_sorted), "species ×", ncol(mat_sorted), "loci"),
    x = "UCE Loci", y = "Species"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 5, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 5),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# 第八步：保存结果 (修改版 - 缩略图预览)
output_dir <- alignment_dir

# 保存热图
png_file <- file.path(output_dir, "completeness_heatmap_preview.png")
pdf_file <- file.path(output_dir, "completeness_heatmap_preview.pdf")

# 设置一个合理的、固定的尺寸（例如 12x8 英寸）
# 这个尺寸在大多数显示器上都能完整显示
FIXED_WIDTH  <- 12
FIXED_HEIGHT <- 8

# 直接保存，不再尝试打印大图
ggsave(png_file, p, 
       width = FIXED_WIDTH, 
       height = FIXED_HEIGHT, 
       dpi = 300, 
       bg = "white",
       device = "png")

ggsave(pdf_file, p, 
       width = FIXED_WIDTH, 
       height = FIXED_HEIGHT, 
       # dpi 在 PDF 中不适用，但可以设置
       bg = "white",
       device = "pdf")