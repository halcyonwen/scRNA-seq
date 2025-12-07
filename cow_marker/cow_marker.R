##########################################################################
### 阶段一：环境准备与数据加载
##########################################################################

# 1. 加载所有必要的程序包
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(monocle3)
library(readxl) 
library(RColorBrewer) # 用于后续病毒绘图配色

# 2. 读取数据
# ------------------------------------------------------------------------
# 对牛乳腺细胞进行分析
seurat_obj <- readRDS("D:/单细胞分析/92501092801_DXF_cowvirus_20250910/MockvsAHvsCXvsTX/MockvsAHvsCXvsTX_combined.rds")
print(seurat_obj)

print("数据加载完成！")
# 绘制初步 UMAP
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
# ------------------------------------------------------------------------
# 牛乳腺细胞 marker 注释
all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(as.data.frame(all_markers), "all_markers.csv")

# 3. 牛乳腺细胞数据过滤
# ------------------------------------------------------------------------
# 过滤掉 3 簇重新绘制 UMAP 图
seurat_obj <- subset(seurat_obj, idents = c("3"), invert = TRUE)
DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("")

##########################################################################
### 阶段三（最终修正版）：自定义顺序 + 多细胞类型强制指定 Marker
##########################################################################

library(dplyr)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

print("开始最终版热图绘制...")

# ========================================================================
# 1. 设定严格的细胞顺序 & 强制显示的基因
# ========================================================================

# 1.1 设定细胞在热图上的排列顺序 (X轴)
my_levels <- c("LumSec", "LumHR", "Basal", "Myoepithelial", "Endothelial", "vSMC", "Fibroblasts")

# 1.2 【新增】设定需要强制显示的基因 (配置列表)
# 格式： "细胞类型名称" = c("基因1", "基因2"...)
# 逻辑： 如果在差异基因里找到了这些基因，优先排在最前；不足10个的位置由其他高表达基因补齐
manual_marker_config <- list(
  "Endothelial" = c("PECAM1", "MOBP", "RAPGEF4"),
  "LumHR"       = c("PGR"),  # 强制添加孕酮受体
  "Fibroblasts" = c("COL6A1", "COL6A2", "COL1A2") # 强制添加胶原蛋白家族
)

# ========================================================================
# 2. 确保 Seurat 对象遵循此顺序
# ========================================================================

# 过滤掉数据中可能不存在的类型
existing_levels <- intersect(my_levels, unique(seurat_obj$major_cell_type))

# 设置因子水平
seurat_obj$major_cell_type <- factor(seurat_obj$major_cell_type, levels = existing_levels)
Idents(seurat_obj) <- "major_cell_type"

# ========================================================================
# 3. 计算差异基因 (如果已计算过则跳过)
# ========================================================================
if (!exists("markers_7_types")) {
  print("正在计算差异基因...")
  markers_7_types <- FindAllMarkers(
    seurat_obj, 
    only.pos = TRUE, 
    min.pct = 0.25, 
    logfc.threshold = 0.25,
    max.cells.per.ident = 500
  )
}

# ========================================================================
# 4. 提取 Top 10 基因 (核心逻辑升级：支持批量强制指定)
# ========================================================================

final_gene_list <- c() # 用于存储最终的基因名顺序
final_df_list <- list() # 用于存储数据框以便检查

# 这里的循环确保了最终基因列表严格按照 existing_levels 的细胞顺序排列
for (ctype in existing_levels) {
  
  # A. 获取当前细胞类型的所有差异基因，按 LogFC 排序
  current_markers <- markers_7_types %>% 
    filter(cluster == ctype) %>% 
    arrange(desc(avg_log2FC))
  
  # B. 检查是否有强制指定的基因
  if (ctype %in% names(manual_marker_config)) {
    # 获取强制基因列表
    force_genes <- manual_marker_config[[ctype]]
    
    # 1. 提取强制基因 (注意：必须是计算出来的差异基因里有的)
    # 如果PGR表达量太低没被FindAllMarkers算出来，这里会自动忽略，防止报错
    df_force <- current_markers %>% filter(gene %in% force_genes)
    
    # 2. 提取剩余的自动基因来补位 (补齐到10个)
    n_needed <- 10 - nrow(df_force)
    df_auto <- current_markers %>% 
      filter(!gene %in% force_genes) %>% 
      head(n = n_needed)
    
    # 3. 合并
    df_final <- bind_rows(df_force, df_auto)
    
  } else {
    # C. 如果没有强制要求，直接取 Top 10
    df_final <- current_markers %>% head(n = 10)
  }
  
  # D. 存储结果
  final_df_list[[ctype]] <- df_final
  final_gene_list <- c(final_gene_list, df_final$gene)
}

# 合并为一个总表 (仅用于查看，不用于画图顺序，画图顺序由 final_gene_list 决定)
final_df_combined <- bind_rows(final_df_list)

print("检查特定基因是否成功入选：")
print(final_df_combined %>% filter(gene %in% c("PGR", "COL6A1", "PECAM1")))

# ========================================================================
# 5. 准备热图数据 (AverageExpression)
# ========================================================================

# 计算平均表达量
avg_exp <- AverageExpression(seurat_obj, group.by = "major_cell_type", layer = "data") 
avg_matrix <- avg_exp$RNA

# 数据清洗与排序
# 1. 筛选基因 (取交集，防止某些强制基因没在矩阵里)
genes_to_keep <- intersect(final_gene_list, rownames(avg_matrix))

# 2. 【关键】按 final_gene_list 的生成顺序排列行 (Y轴顺序)
# 这里必须用 match 确保顺序完全一致
avg_matrix <- avg_matrix[match(genes_to_keep, rownames(avg_matrix)), ] 

# 3. 【关键】按 existing_levels 的顺序排列列 (X轴顺序)
cols_to_keep <- intersect(existing_levels, colnames(avg_matrix))
avg_matrix <- avg_matrix[, cols_to_keep] 

# Z-score 标准化
plot_matrix <- t(scale(t(avg_matrix)))
plot_matrix[is.nan(plot_matrix)] <- 0
plot_matrix[plot_matrix > 2] <- 2
plot_matrix[plot_matrix < -2] <- -2

# ========================================================================
# 6. 绘制并保存热图
# ========================================================================

# 配色方案
col_fun <- colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026"))

# 保存为 PDF
pdf("Average_Heatmap_Custom_Markers.pdf", width = 6, height = 10)

ht <- Heatmap(plot_matrix,
              name = "Z-Score",          
              col = col_fun,            
              
              # --- 布局设置 ---
              row_names_side = "left",   
              column_names_side = "bottom", 
              column_names_rot = 45,       
              
              # --- 排序控制 (全部关闭，依赖数据本身的顺序) ---
              cluster_rows = FALSE,      
              cluster_columns = FALSE,   
              
              # --- 美化细节 ---
              rect_gp = gpar(col = "white", lwd = 1), 
              row_names_gp = gpar(fontsize = 10),      
              column_names_gp = gpar(fontsize = 12, fontface = "bold"), 
              
              # --- 增加分割线 ---
              column_split = factor(colnames(plot_matrix), levels = existing_levels),
              column_gap = unit(0.2, "mm"),
              border = TRUE
)

draw(ht)
dev.off()

print("绘图完成！文件已保存为 Average_Heatmap_Custom_Markers.pdf")
