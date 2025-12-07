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
### 阶段三：牛乳腺细胞核心注释与热图 (Heatmap)
##########################################################################

# 1. 大类划分逻辑 (Cluster -> Cell Type)
# ------------------------------------------------------------------------
cluster_to_celltype_map <- list(
  "LumSec" = c( 0,2,14 ),
  "LumHR" = c( 15 ),
  "Fibroblasts" = c( 7,13,16 ), 
  "Basal" = c( 1,4,6,10,11,12 ),
  "Endothelial" = c(5),
  "vSMC" = c(8),
  "Myoepithelial" = c(9)
)

# 检查依赖列
if(!"seurat_clusters" %in% colnames(seurat_obj@meta.data)){
  stop("错误：找不到 'seurat_clusters' 列，请确保已完成聚类分析(FindClusters)。")
}

# 创建新的元数据列并赋值
seurat_obj$major_cell_type <- as.character(seurat_obj$seurat_clusters)

for (cell_type in names(cluster_to_celltype_map)) {
  cluster_ids <- cluster_to_celltype_map[[cell_type]]
  cells_to_rename <- seurat_obj$seurat_clusters %in% cluster_ids
  seurat_obj$major_cell_type[cells_to_rename] <- cell_type
}

# 设置身份与顺序
Idents(seurat_obj) <- "major_cell_type"
seurat_obj$major_cell_type <- factor(
  seurat_obj$major_cell_type, 
  levels = names(cluster_to_celltype_map)
)
Idents(seurat_obj) <- seurat_obj$major_cell_type

# 检查并可视化
print(table(seurat_obj$major_cell_type))
DimPlot(seurat_obj, reduction = "umap", group.by = "major_cell_type", label = TRUE) +
  ggtitle("Major Cell Types of Bovine Mammary Gland")

##########################################################################
### 阶段三（最终修正版）：自定义顺序 + 强制内皮基因
##########################################################################

library(dplyr)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

print("开始最终版热图绘制...")

# 1. 设定严格的细胞顺序 (您的核心需求)
# ------------------------------------------------------------------------
my_levels <- c("LumSec", "LumHR", "Basal", "Myoepithelial", "Endothelial", "vSMC", "Fibroblasts")

# 2. 确保 Seurat 对象遵循此顺序
# ------------------------------------------------------------------------
# 过滤掉数据中可能不存在的类型
existing_levels <- intersect(my_levels, unique(seurat_obj$major_cell_type))

# 设置因子水平 (这将决定热图列的顺序)
seurat_obj$major_cell_type <- factor(seurat_obj$major_cell_type, levels = existing_levels)
Idents(seurat_obj) <- "major_cell_type"

# 3. 计算差异基因 (如果已计算过 markers_7_types，可跳过此步)
# ------------------------------------------------------------------------
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
# 4. 提取 Top 10 基因 (核心逻辑：强制内皮基因 + 顺序对齐)
# ========================================================================

# A. 定义内皮细胞必须包含的基因
endo_force_genes <- c("PECAM1", "MOBP", "RAPGEF4")

# B. 拆分处理

# --- Part 1: 非内皮细胞 (标准 Top 10) ---
top10_others <- markers_7_types %>%
  filter(cluster != "Endothelial") %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  ungroup() 

# --- Part 2: 内皮细胞 (强制包含特定基因) ---
endo_all_markers <- markers_7_types %>% 
  filter(cluster == "Endothelial") %>%
  arrange(desc(avg_log2FC))

# (1) 提取强制基因
valid_force <- endo_all_markers %>% 
  filter(gene %in% endo_force_genes)

# (2) 提取自动基因 (补齐到10个)
auto_genes <- endo_all_markers %>% 
  filter(!gene %in% endo_force_genes) %>%
  head(n = 10 - nrow(valid_force)) 

# (3) 合并内皮基因 (强制基因在前)
top10_endo <- bind_rows(valid_force, auto_genes)

# C. 合并所有数据
final_df <- bind_rows(top10_others, top10_endo)

# D. 【关键步骤】对最终列表按 my_levels 重新排序
# 如果不加这一步，Endothelial的基因可能会跑到列表最后，导致热图错位
final_df$cluster <- factor(final_df$cluster, levels = existing_levels)
final_df <- final_df %>% arrange(cluster) 

# 提取最终基因名
final_gene_list <- final_df$gene

print("检查内皮细胞基因顺序 (应包含 PECAM1):")
print(final_gene_list[final_df$cluster == "Endothelial"])

# ========================================================================
# 5. 准备热图数据 (AverageExpression)
# ========================================================================

# 计算平均表达量
avg_exp <- AverageExpression(seurat_obj, group.by = "major_cell_type", layer = "data") 
avg_matrix <- avg_exp$RNA

# 数据清洗与排序 (确保行和列都严格对齐)
# 1. 筛选基因
genes_to_keep <- intersect(final_gene_list, rownames(avg_matrix))
avg_matrix <- avg_matrix[genes_to_keep, ]

# 2. 【关键】按 final_gene_list 的顺序排列行 (Y轴顺序)
avg_matrix <- avg_matrix[match(genes_to_keep, rownames(avg_matrix)), ] 

# 3. 【关键】按 my_levels 的顺序排列列 (X轴顺序)
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

# 配色方案 (经典红白蓝)
col_fun <- colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026"))

# 保存为 PDF
pdf("Average_Heatmap_Fixed_Order.pdf", width = 6, height = 10)

ht <- Heatmap(plot_matrix,
              name = "Z-Score",          
              col = col_fun,            
              
              # --- 布局设置 ---
              row_names_side = "left",   
              column_names_side = "bottom", 
              column_names_rot = 45,      
              
              # --- 排序控制 (全部关闭，完全听从我们的指令) ---
              cluster_rows = FALSE,      
              cluster_columns = FALSE,   
              
              # --- 美化细节 ---
              rect_gp = gpar(col = "white", lwd = 1), 
              row_names_gp = gpar(fontsize = 10),      
              column_names_gp = gpar(fontsize = 12, fontface = "bold"), 
              
              # --- 增加分割线 (可选) ---
              # 这会在不同细胞类型之间画一条粗线，让区分更明显
              column_split = factor(colnames(plot_matrix), levels = existing_levels),
              column_gap = unit(0.2, "mm"),
              border = TRUE
)

draw(ht)
dev.off()

print("绘图完成！文件已保存为 Average_Heatmap_Fixed_Order.pdf")
