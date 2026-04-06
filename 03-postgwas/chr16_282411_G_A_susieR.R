# 从CRAN安装
install.packages("susieR")

# 或从GitHub安装开发版
# devtools::install_github("stephenslab/susieR")

# 加载包
rm(list=ls())
library(susieR)
library(data.table)
eq=read.table("C://new_H_disease_NatureGenetics//Fig4-postGWAS//finemap-SUSIER//case-control//chr16_282411_G_A.zscore.txt", sep="\t", stringsAsFactors=FALSE,header=T)
dat <- fread("C://new_H_disease_NatureGenetics//Fig4-postGWAS//finemap-SUSIER//case-control//chr16.LD.matrix.ld")
head(eq)
dat1=as.matrix(dat)
fitted_rss3 <- susie_rss(c(eq$zscore), dat1, n=1831, L = 5)
susie_plot(fitted_rss3, y="PIP", add_legend=TRUE,
           cex.axis = 1.5,    # 坐标轴字体增大50%
           cex.lab = 1.5,     # 标签字体增大50%
           cex = 1.5)      # 普通点大小（默认1.0）
          
summary(fitted_rss3)$cs
top_snps <- c(160, 217, 214, 276, 102, 748, 516, 564, 136, 265, 268)
# 确认这些SNP的PIP是否接近1
top_pips <- fitted_rss3$pip[top_snps]
data.frame(SNP=top_snps, PIP=top_pips)
eq[top_snps, ]  

