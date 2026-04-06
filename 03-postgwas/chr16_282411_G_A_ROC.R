# 设置镜像并安装
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
install.packages("survminer", dependencies = TRUE)
# 强制安装最底层的缺失包
install.packages("gridtext")
# 顺便把可能也缺的 ggtext 一起重装/补齐
install.packages("ggtext")
# ==========================================
# 1. 环境准备与数据加载
# ==========================================
rm(list = ls())
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

# 读取数据 (建议使用 check.names=FALSE 保持基因型列名原始格式)
file_path <- "C://new_H_disease_NatureGenetics//Fig4-postGWAS//genotype-phenotype//case-control//chr16_282411_G_A_ROC_Survival time without transfusion.txt"
data2 <- read.table(file_path, header = TRUE, sep = "\t", check.names = FALSE)

# ==========================================
# 2. 数据清洗与变量预处理
# ==========================================
data2 <- data2 %>%
  # 确保数值型转换
  mutate(
    Time = as.numeric(Survivaltimewithouttransfusionmonths),
    Status = as.numeric(Status)
  ) %>%
  # 剔除关键变量缺失的样本
  filter(!is.na(Time), !is.na(Status), !is.na(chr16_282411_G_A)) %>%
  # 【关键】将 GG 设置为参考组 (Reference)，以便计算相对于 GG 的风险
  mutate(
    chr16_282411_G_A = factor(chr16_282411_G_A, levels = c("GG", "GA", "AA")),
    # 创建合并分组：GG vs. Non-GG (GA+AA)
    Group_Combined = if_else(chr16_282411_G_A == "GG", "GG (Protected)", "Non-GG (Risk)"),
    Group_Combined = factor(Group_Combined, levels = c("GG (Protected)", "Non-GG (Risk)"))
  )

# 快速检查状态分布
print(table(data2$Status, data2$chr16_282411_G_A))

# ==========================================
# 3. 统计建模：Cox 比例风险模型 (Hazard Ratio)
# ==========================================
# 计算原始三组的 HR
cox_fit <- coxph(Surv(Time, Status) ~ chr16_282411_G_A, data = data2)
summary(cox_fit) # 重点关注 exp(coef) 和对应的 95% CI

# 计算合并组的 HR (Non-GG vs GG)
cox_combined <- coxph(Surv(Time, Status) ~ Group_Combined, data = data2)
summary(cox_combined)

# ==========================================
# 4. 高级可视化：GG vs Non-GG (JCO 风格)
# ==========================================
fit_combined <- survfit(Surv(Time, Status) ~ Group_Combined, data = data2)

# 5. 执行美化后的绘图
p_final <- ggsurvplot(
  fit_combined, 
  data = data2,
  
  # --- 1. 统计方法与P值设置 ---
  pval = FALSE,                # 不显示 P 值
  pval.method = TRUE,         # 显示 "Log-rank" 字样
  pval.coord = c(15, 0.1),   # 【关键】修改P值坐标(x, y)，将其移到左下角空白处，避免重叠
  pval.size = 6,              # 加大P值字体
  
  # --- 2. 曲线与置信区间 ---
  conf.int = TRUE,            # 显示置信区间
  palette = "npg",            # JCO 风格配色
  size = 1.2,                 # 增加线条粗细
  
  # --- 3. 标签与标题 ---
  legend.title = "Genotype",
  legend.labs = c("GG (Ref)", "GA/AA (Risk)"),
  xlab = "Follow-up Time (Months)", 
  ylab = "Transfusion-Free Probability",
  title = "chr16:282411_G/A",
  # --- 4. 刻度与风险表 ---
  break.time.by = 60,         # 每5年一个刻度
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.y.text = FALSE,  # 风险表左侧不显示文字，更整洁
  risk.table.height = 0.25,
  
  # --- 5. 样式深度定制 (字体加大/背景空白) ---
  ggtheme = theme_classic() 
)
p_final$plot <- p_final$plot + 
  annotate(
    "text", 
    x = 30, y = 0.10,          # 【在这里调坐标】x 是横轴月份，y 是纵轴概率
    label = "Log-rank\np = 0.00036", 
    size = 5,                 # 字体大小
    hjust = 0,                # 左对齐
    fontface = "bold",        # 加粗
    lineheight = 1.1          # 控制 Log-rank 和 p 之间的行间距
  ) +
  theme(
    axis.line = element_line(linewidth = 1.2, color = "black"),
    axis.ticks = element_line(linewidth = 1.2, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
  )

# 修改风险表的字体（ggsurvplot 的风险表需要单独微调）
p_final$table <- p_final$table + 
  theme(
    axis.line = element_line(linewidth = 1.0, color = "black"),
    axis.ticks.x = element_line(linewidth = 1.0, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold")
  )
# --- 输出预览 ---
print(p_final)

# 导出高清投稿图 (建议 PDF，300 DPI 效果最佳)
pdf("C://new_H_disease_NatureGenetics//Fig4-postGWAS//genotype-phenotype//case-control//chr16_282411_G_A_ROC-1.pdf", 
    width = 8, height = 8, onefile = FALSE)
# 2. 打印你的图片对象
print(p_final)
# 3. 关闭画布（这步极其重要，否则 PDF 文件会损坏无法打开）
dev.off()

