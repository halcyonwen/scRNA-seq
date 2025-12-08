# 步骤0: 加载必要的R包
library(Seurat)
library(ggplot2)
library(dplyr)

# 假设您的Seurat对象是 sc_object
sc_object <- readRDS("D:/单细胞分析/92501092801_DXF_cowvirus_20250910/MockvsAHvsCXvsTX/MockvsAHvsCXvsTX_combined.rds")

# --- 步骤1: 定义“病毒表达量” ---
# 我们将优先使用方法B（计算特征分数），因为它更稳健。

# 首先，定义一个包含所有病毒基因名称的列表（vector）
# !!! 关键：您需要根据您的数据集，将这里的基因名替换为真实的病毒基因名 !!!
viral_genes <- c("CX-virus-M") 

# 检查这些基因是否存在于您的数据中
viral_genes_in_data <- viral_genes[viral_genes %in% rownames(sc_object)]
if (length(viral_genes_in_data) == 0) {
  stop("错误：您定义的病毒基因在数据集中一个也找不到，请检查基因名称。")
}

# 使用 AddModuleScore 计算病毒特征分数
# a new column named 'ViralScore1' will be added to sc_object@meta.data
sc_object <- AddModuleScore(
  object = sc_object,
  features = list(viral_genes_in_data),
  name = 'ViralScore' # 生成的metadata列名是 ViralScore1
)

# --- 步骤2 & 3: 提取数据并整合筛选 ---

# 从Seurat对象中提取需要的信息
# 1. Mx1 基因的表达数据
mx1_expression <- GetAssayData(sc_object, assay = "RNA", slot = "data")["MX1", ]

# 2. 细胞分组信息 和 我们刚计算的病毒分数
metadata <- sc_object@meta.data

# --- 步骤 1: 创建原始数据框 ---

plot_df <- data.frame(
  Mx1 = mx1_expression,
  Viral_Score = sc_object@meta.data[["ViralScore1"]],
  Group = sc_object@meta.data[["stim"]]
)

# --- 步骤 2: 创建新的、合并后的分组列 ---
# 我们将添加一个名为 'condition' 的新列
# gsub("\\d+$", "", Group) 的作用是找到(Group)列中每个字符串末尾($)的一个或多个数字(\\d+)并替换为空("")
# 例如，"AH1" 会变成 "AH"，"CX2" 会变成 "CX"
library(dplyr)
plot_df <- plot_df %>%
  mutate(condition = gsub("\\d+$", "", Group))

# 检查一下新列是否正确创建
print(head(plot_df))


# --- 步骤 3: 使用新的'condition'列进行筛选和绘图 ---
# 定义我们感兴趣的条件
target_conditions <- c("AH", "CX", "TX")

# 使用新的 'condition' 列进行筛选
plot_df_filtered <- plot_df %>%
  filter(condition %in% target_conditions)

# 确认筛选成功
print(paste("成功筛选数据，找到", nrow(plot_df_filtered), "个细胞用于绘图。"))


# --- 步骤 4: 可视化 ---
# 注意，现在 facet_wrap 使用的是新的 '~ condition' 列
correlation_plot <- ggplot(plot_df_filtered, aes(x = Viral_Score, y = Mx1)) +
  geom_point(alpha = 0.5, color = "orange") +
  geom_smooth(method = "lm", se = FALSE, color = "purple") +
  facet_wrap(~ condition, scales = "free") + # <-- 使用新的 condition 列分面
  labs(
    title = "Mx1 Expression vs. Viral Score",
    x = "Viral Expression Score",
    y = "Mx1 Expression Level"
  ) +
  theme_classic()

# 显示图形
print(correlation_plot)

# --- 重新运行从步骤4开始的所有代码以确保所有对象都存在 ---

# 假设 plot_df 已经成功创建并包含了 'Mx1', 'Viral_Score', 'Group', 'condition' 列

# --- 步骤 5: 筛选数据 ---
target_conditions <- c("AH", "CX", "TX")
plot_df_filtered <- plot_df %>%
  filter(condition %in% target_conditions)

# --- 步骤 6: 统计学计算 (新步骤) ---
# 我们将按 'condition' 分组，对每组计算相关性，并创建一个用于标注的新数据框

library(dplyr) # 确保dplyr已加载


stats_summary <- plot_df_filtered %>%
  group_by(condition) %>%
  summarise(
    # 我们仍然运行一次cor.test来方便地获取rho值
    rho = cor(Viral_Score, Mx1, method = "spearman"),
    # 同时，我们计算每个分组的细胞数n
    n = n(),
    .groups = 'drop'
  ) %>%
  # --- 手动计算t统计量 ---
  mutate(
    # 根据公式计算t值
    t_statistic = rho * sqrt((n - 2) / (1 - rho^2)),
    # 自由度是 n - 2
    df = n - 2
  ) %>%
  # --- 高精度P值计算 ---
  mutate(
    # 现在使用我们自己计算的 t_statistic 和 df
    log_p_value = log(2) + pt(abs(t_statistic), df = df, lower.tail = FALSE, log.p = TRUE),
    
    # 转换为以10为底的对数
    log10_p_value = log_p_value / log(10),
    
    # 分离指数和尾数
    exponent = floor(log10_p_value),
    mantissa = 10^(log10_p_value - exponent),
    
    # 创建最终的文本标签
    label_text = sprintf("ρ = %.2f\nP ≈ %.2f x 10^%d", rho, mantissa, exponent)
  )

# 查看最终的统计结果数据框
print(stats_summary)


# --- 步骤 7: 可视化 ---
correlation_plot_annotated <- ggplot(plot_df_filtered, aes(x = Viral_Score, y = Mx1)) +
  geom_point(alpha = 0.5, color = "orange", shape = 16, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "purple") +
  facet_wrap(~ condition, scales = "free") +
  
  geom_text(
    data = stats_summary,
    aes(x = Inf, y = Inf, label = label_text),
    hjust = 1.1,
    vjust = 1.5,
    size = 4,
    inherit.aes = FALSE
  ) +
  
  labs(
    title = "Mx1 Expression vs. Viral Score by Condition",
    x = "Viral Expression Score (Unitless Module Score)",
    y = "Mx1 Expression (Log-Normalized)"
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# 显示最终的、标注正确的图形
print(correlation_plot_annotated)
###############################################################################
############################UMAP图展示Mx1基因在不同细胞簇的表达情况############
# 使用 split.by 按您的分组（例如 stim）拆分，查看差异
# 这将清楚地显示Mx1是否在特定刺激组中上调
FeaturePlot(sc_object, 
            features = "MX1", 
            split.by = "stim", # 使用您的分组列
            cols = c("grey", "red") # 定义颜色
)
###############################################################################
###############################################################################
