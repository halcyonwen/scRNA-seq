##########################################################################
### 阶段一：环境准备与数据加载
##########################################################################

# 1. 加载所有必要的程序包
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(readxl) 
library(ComplexHeatmap) 
library(circlize)
library(RColorBrewer)
library(scales) # 用于颜色缩放

# 2. 读取数据 (人乳腺细胞)
# ------------------------------------------------------------------------
seurat_obj <- readRDS("D:/单细胞分析/92501092801_DXF_homo_B1_20250910/MockvsAHvsCXvsTX/MockvsAHvsCXvsTX_combined.rds")
print(seurat_obj)

# ========================================================================
# 1. 再次确认基因列表存在性 (防止 FetchData 报错)

# --- 步骤 1: 识别病毒基因 ---
all_genes <- rownames(seurat_obj)
viral_patterns <- c("^CX-", "^AH-", "^TX-") 
combined_pattern <- paste(viral_patterns, collapse = "|")
viral_genes_found <- grep(combined_pattern, all_genes, value = TRUE, ignore.case = TRUE)

# 只有在找到了病毒基因时才进行后续操作
if (length(viral_genes_found) > 0) {
  
  #################################################################
  # --- 步骤 2: 简化的自定义基因排序逻辑 ---
  #################################################################
  
  # 1. 定义基因排序规则
  major_order <- c("^CX-", "^AH-", "^TX-")
  internal_order <- c("PB2", "PB1", "PA", "NP", "HA", "NA", "M", "NS")
  
  # 2. 构造搜索模式
  all_patterns_matrix <- outer(major_order, internal_order, function(major, internal) {
    paste0(major, ".*", internal)
  })
  ordered_patterns <- as.vector(t(all_patterns_matrix))
  
  # 3. 匹配并组合
  ordered_genes_list <- lapply(ordered_patterns, function(pattern) {
    grep(pattern, viral_genes_found, value = TRUE, ignore.case = TRUE)
  })
  ordered_viral_genes <- unique(unlist(ordered_genes_list))
  
  # 4. 处理剩余基因
  remaining_genes <- setdiff(viral_genes_found, ordered_viral_genes)
  final_ordered_genes <- c(ordered_viral_genes, remaining_genes)
}
# ========================================================================
# 确保所有要画的基因都在对象里
valid_genes <- intersect(final_ordered_genes, rownames(seurat_obj))

if (length(valid_genes) > 0) {
  
  # ========================================================================
  # 2. 手动提取数据 (Bypass Seurat DotPlot)
  # ========================================================================
  print("正在提取表达矩阵数据...")
  
  # 提取基因表达数据和细胞分类
  # FetchData 返回的是一个巨大的表格：行是细胞，列是基因 + ident
  plot_data <- FetchData(seurat_obj, vars = c(valid_genes, "ident"))
  
  # 转换为长格式 (Tidy format)，方便 ggplot 处理
  # 这一步如果不报错，说明数据就是干净的
  plot_data_long <- plot_data %>%
    as.data.frame() %>%
    pivot_longer(
      cols = -ident,        # 除了 ident 列，其他都是基因列
      names_to = "Gene",    # 列名变为 Gene 字段
      values_to = "Expression" # 数值变为 Expression 字段
    )
} else {
  # ========================================================================
  # 3. 计算统计量 (复刻 DotPlot 的计算逻辑)
  # ========================================================================
  print("正在计算平均表达量和表达比例...")
  
  # DotPlot 的两个核心指标：
  # 1. Pct: 该群中有多少比例的细胞表达量 > 0
  # 2. Avg: 该群中细胞的平均表达量
  
  summary_data <- plot_data_long %>%
    group_by(ident, Gene) %>%
    summarise(
      Pct = sum(Expression > 0) / n() * 100,  # 表达比例 (%)
      Avg = mean(expm1(Expression)),          # 平均表达量 (反log处理，Seurat默认存的是log值)
      .groups = 'drop'
    )
  cell_type_list <- list(
    "Myeloid"        = c("0", "6", "7","8","9","10","11"),
    "Fibroblast"         = c("1","2","5","12","13","14","15","16","17","18"),
    "Basal"         = c("3","4"),
    "LumSec" = c("19")
  )
  
  new_cluster_order <- unlist(cell_type_list)
  
  # --- 应用顺序 ---
  valid_levels <- intersect(new_cluster_order, unique(as.character(Idents(seurat_obj))))
  seurat_obj@active.ident <- factor(seurat_obj@active.ident, levels = valid_levels)
  # ========================================================================
  # 4. 数据标准化 (Z-score) - 为了画出好看的红蓝色
  # ========================================================================
  # 对每个基因在不同 Cluster 之间进行 Z-score 标准化
  summary_data <- summary_data %>%
    group_by(Gene) %>%
    mutate(
      ScaleAvg = scale(Avg), # 计算 Z-score
      # 截断极值，防止颜色过爆 (类似于 Seurat 的 min/max cutoff)
      ScaleAvg = pmin(pmax(ScaleAvg, -2.5), 2.5) 
    ) %>%
    ungroup()
  # ========================================================================
  # 4. 动态计算分割线位置
  # ========================================================================
  separator_positions <- c()
  current_pos <- 0
  
  for (i in 1:(length(cell_type_list) - 1)) { 
    group_clusters <- cell_type_list[[i]]
    n_present <- sum(group_clusters %in% valid_levels)
    if (n_present > 0) {
      current_pos <- current_pos + n_present
      separator_positions <- c(separator_positions, current_pos + 0.5)
    }
  }
  # ========================================================================
  # 5. 强制设定顺序 (Factor Reordering)
  # ========================================================================
  # 基因顺序 (X轴)
  summary_data$Gene <- factor(summary_data$Gene, levels = valid_genes)
  
  # 细胞群顺序 (Y轴) - 之前你已经设定好了 seurat_obj 的 levels，这里直接继承
  # 注意：DotPlot里用了 coord_flip，所以这里我们要小心对应 X 和 Y
  summary_data$ident <- factor(summary_data$ident, levels = levels(seurat_obj))
  
  # ========================================================================
  # 6. 使用 ggplot2 绘图 (完全复刻 DotPlot 样式)
  # ========================================================================
  print("正在绘制最终图表...")
  
  # 设定颜色：低=蓝，中=白，高=红
  my_gradient <- c("#313695", "white", "#A50026")
  
  p_manual_dot <- ggplot(summary_data, aes(x = Gene, y = ident)) +
    # 画点：大小代表比例，颜色代表相对表达量
    geom_point(aes(size = Pct, color = ScaleAvg)) + 
    
    # 颜色标尺
    scale_color_gradientn(colors = my_gradient, name = "Avg Exp\n(Z-score)") +
    
    # 大小标尺
    scale_size(range = c(0, 6), name = "Percent\nExpressed") +
    
    # 翻转坐标轴 (X轴变Cluster，Y轴变Gene，符合你之前 DotPlot + coord_flip 的效果)
    coord_flip() + 
    
    # 【关键】添加虚线分割 (使用你之前算好的 separator_positions)
    geom_hline(yintercept = separator_positions, 
               linetype = "dashed", color = "black", linewidth = 0.8) +
    
    # 主题美化
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(face = "italic", color = "black"),
      panel.grid.major = element_line(color = "grey95"), # 淡网格线
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
    ) +
    
    labs(x = "Viral Genes", y = "Cell Type Clusters", 
         title = "Viral Gene Expression (Manually Constructed)")
  
  # ========================================================================
  # 7. 保存结果
  # ========================================================================
  calc_height <- 2 + (length(valid_genes) * 0.35)
  
  ggsave("viral_expression_dotplot_human_MANUAL.pdf", p_manual_dot, width = 14, height = calc_height, limitsize = FALSE)
  
  print(">>> 成功！已生成: viral_expression_dotplot_cow_MANUAL.pdf")
  print(p_manual_dot)
}  else {
    print("错误：经过过滤后没有有效的病毒基因可供绘图。")
  }
