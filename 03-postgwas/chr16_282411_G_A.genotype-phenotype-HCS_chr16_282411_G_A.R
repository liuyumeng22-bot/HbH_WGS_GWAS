install.packages("patchwork")
install.packages("ggsci")
install.packages("readxl")

rm(list = ls())
library(ggpubr) 
library(patchwork) 
library(ggsci)
library(tidyverse)
library(readxl)
library(ggplot2)
data1=read_table("C://new_H_disease_NatureGenetics//Fig4-postGWAS//genotype-phenotype//case-control//chr16_282411_G_A.genotype-phenotype-HCS.txt")
head(data1)
str(data1)
data2 <- mutate(data1, across(3:9, as.numeric))

compare_list1 <- list(
  c("AA","GG"),c("GA","GG"),c("AA","GA"))

p1 <- ggplot(data= data2,aes(x = factor(chr16_282411_G_A,levels =c("GG","GA","AA")),y =Hemoglobin)) +
  geom_boxplot(aes(fill = chr16_282411_G_A),alpha = 0.6,outlier.color=NA)+ 
  geom_jitter(width=0.2,alpha = 0.3)+
  scale_fill_manual(values = c("#f79f1f","#a3cb38","#1289a7"))+
  stat_compare_means(
    comparisons = compare_list1,
    method = "wilcox.test",
    label = "p.signif",
    size=6)+
  theme_bw()+
  theme(legend.position="none", panel.grid = element_blank(),
        axis.text.x=element_text(colour="black",size=16), #????x???潭缺?签??????????
        axis.text.y=element_text(size=16), #????y???潭缺?签??????????
        axis.title.y=element_text(size = 18), #????y???谋?????????????
        axis.title.x=element_text(size = 18),)+ 
  scale_y_continuous(expand = c(0.1, 0.1))+
  ylab ("Hemoglobin(g/L)")+
  xlab("chr16:282411_G/A)\n(rs199594877)")

p1
ggsave(plot = p1, filename = "C://new_H_disease_NatureGenetics//Fig4-postGWAS//genotype-phenotype//case-control//HCS-chr16_282411_G_A_Hemoglobin_pheno_snp_assoc.pdf", width = 6, height = 6)

################################
compare_list1 <- list(
  c("AA","GG"),c("GA","GG"),c("AA","GA"))

p1 <- ggplot(data= data2,aes(x = factor(chr16_282411_G_A,levels =c("GG","GA","AA")),y =TBIL)) +
  geom_boxplot(aes(fill = chr16_282411_G_A),alpha = 0.6,outlier.color=NA)+ 
  geom_jitter(width=0.2,alpha = 0.3)+
  scale_fill_manual(values = c("#f79f1f","#a3cb38","#1289a7"))+
  stat_compare_means(
    comparisons = compare_list1,
    method = "wilcox.test",
    label = "p.signif",
    size=6,
    method.args = list(exact = FALSE))+
  theme_bw()+
  theme(legend.position="none", panel.grid = element_blank(),
        axis.text.x=element_text(colour="black",size=16), #????x???潭缺?签??????????
        axis.text.y=element_text(size=16), #????y???潭缺?签??????????
        axis.title.y=element_text(size = 18), #????y???谋?????????????
        axis.title.x=element_text(size = 18),)+ 
  scale_y_continuous(expand = c(0.1, 0.1))+
  ylab ("Total bilirubin(μmol/L)")+
  xlab("chr16:282411_G/A)\n(rs199594877)")

p1
ggsave(plot = p1, filename = "C://new_H_disease_NatureGenetics//Fig4-postGWAS//genotype-phenotype//case-control//HCS-chr16_282411_G_A_TBIL_pheno_snp_assoc.pdf", width = 6, height = 6)


################################
compare_list1 <- list(
  c("AA","GG"),c("GA","GG"),c("AA","GA"))

p1 <- ggplot(data= data2,aes(x = factor(chr16_282411_G_A,levels =c("GG","GA","AA")),y =RBC)) +
  geom_boxplot(aes(fill = chr16_282411_G_A),alpha = 0.6,outlier.color=NA)+ 
  geom_jitter(width=0.2,alpha = 0.3)+
  scale_fill_manual(values = c("#f79f1f","#a3cb38","#1289a7"))+
  stat_compare_means(
    comparisons = compare_list1,
    method = "wilcox.test",
    label = "p.signif",
    size=6,
    method.args = list(exact = FALSE))+
  theme_bw()+
  theme(legend.position="none", panel.grid = element_blank(),
        axis.text.x=element_text(colour="black",size=16), #????x???潭缺?签??????????
        axis.text.y=element_text(size=16), #????y???潭缺?签??????????
        axis.title.y=element_text(size = 18), #????y???谋?????????????
        axis.title.x=element_text(size = 18),)+ 
  scale_y_continuous(expand = c(0.1, 0.1))+
  ylab ("Red blood cell count(E12/L)")+
  xlab("chr16:282411_G/A)\n(rs199594877)")

p1
ggsave(plot = p1, filename = "C://new_H_disease_NatureGenetics//Fig4-postGWAS//genotype-phenotype//case-control//HCS-chr16_282411_G_A_RBC_pheno_snp_assoc.pdf", width = 6, height = 6)

################################
compare_list1 <- list(
  c("AA","GG"),c("GA","GG"),c("AA","GA"))

p1 <- ggplot(data= data2,aes(x = factor(chr16_282411_G_A,levels =c("GG","GA","AA")),y =Hematocrit)) +
  geom_boxplot(aes(fill = chr16_282411_G_A),alpha = 0.6,outlier.color=NA)+ 
  geom_jitter(width=0.2,alpha = 0.3)+
  scale_fill_manual(values = c("#f79f1f","#a3cb38","#1289a7"))+
  stat_compare_means(
    comparisons = compare_list1,
    method = "wilcox.test",
    label = "p.signif",
    size=6,
    method.args = list(exact = FALSE))+
  theme_bw()+
  theme(legend.position="none", panel.grid = element_blank(),
        axis.text.x=element_text(colour="black",size=16), #????x???潭缺?签??????????
        axis.text.y=element_text(size=16), #????y???潭缺?签??????????
        axis.title.y=element_text(size = 18), #????y???谋?????????????
        axis.title.x=element_text(size = 18),)+ 
  scale_y_continuous(expand = c(0.1, 0.1))+
  ylab ("Hematocrit (L/L)")+
  xlab("chr16:282411_G/A)\n(rs199594877)")

p1
ggsave(plot = p1, filename = "C://new_H_disease_NatureGenetics//Fig4-postGWAS//genotype-phenotype//case-control//HCS-chr16_282411_G_A_Hematocrit_pheno_snp_assoc.pdf", width = 6, height = 6)

################################
compare_list1 <- list(
  c("AA","GG"),c("GA","GG"),c("AA","GA"))

p1 <- ggplot(data= data2,aes(x = factor(chr16_282411_G_A,levels =c("GG","GA","AA")),y =Reticulocyte)) +
  geom_boxplot(aes(fill = chr16_282411_G_A),alpha = 0.6,outlier.color=NA)+ 
  geom_jitter(width=0.2,alpha = 0.3)+
  scale_fill_manual(values = c("#f79f1f","#a3cb38","#1289a7"))+
  stat_compare_means(
    comparisons = compare_list1,
    method = "wilcox.test",
    label = "p.signif",
    size=6,
    method.args = list(exact = FALSE))+
  theme_bw()+
  theme(legend.position="none", panel.grid = element_blank(),
        axis.text.x=element_text(colour="black",size=16), #????x???潭缺?签??????????
        axis.text.y=element_text(size=16), #????y???潭缺?签??????????
        axis.title.y=element_text(size = 18), #????y???谋?????????????
        axis.title.x=element_text(size = 18),)+ 
  scale_y_continuous(expand = c(0.1, 0.1))+
  ylab ("Reticulocyte count(E12/L)")+
  xlab("chr16:282411_G/A)\n(rs199594877)")

p1
ggsave(plot = p1, filename = "C://new_H_disease_NatureGenetics//Fig4-postGWAS//genotype-phenotype//case-control//HCS-chr16_282411_G_A_Reticulocyte_pheno_snp_assoc.pdf", width = 6, height = 6)

################################
compare_list1 <- list(
  c("AA","GG"),c("GA","GG"),c("AA","GA"))

p1 <- ggplot(data= data2,aes(x = factor(chr16_282411_G_A,levels =c("GG","GA","AA")),y =Ferritin)) +
  geom_boxplot(aes(fill = chr16_282411_G_A),alpha = 0.6,outlier.color=NA)+ 
  geom_jitter(width=0.2,alpha = 0.3)+
  scale_fill_manual(values = c("#f79f1f","#a3cb38","#1289a7"))+
  stat_compare_means(
    comparisons = compare_list1,
    method = "wilcox.test",
    label = "p.signif",
    size=6,
    method.args = list(exact = FALSE))+
  theme_bw()+
  theme(legend.position="none", panel.grid = element_blank(),
        axis.text.x=element_text(colour="black",size=16), #????x???潭缺?签??????????
        axis.text.y=element_text(size=16), #????y???潭缺?签??????????
        axis.title.y=element_text(size = 18), #????y???谋?????????????
        axis.title.x=element_text(size = 18),)+ 
  scale_y_continuous(expand = c(0.1, 0.1))+
  ylab ("Ferritin(ng/mL)")+
  xlab("chr16:282411_G/A)\n(rs199594877)")

p1
ggsave(plot = p1, filename = "C://new_H_disease_NatureGenetics//Fig4-postGWAS//genotype-phenotype//case-control//HCS-chr16_282411_G_A_Ferritin_pheno_snp_assoc.pdf", width = 6, height = 6)


################################
################################

p1 <- ggplot(data= data2,aes(x = factor(chr16_282411_G_A,levels =c("GG","GA","AA")),y =Survivaltimewithouttransfusionmonths)) +
  geom_boxplot(aes(fill = chr16_282411_G_A),alpha = 0.6,outlier.color=NA)+ 
  geom_jitter(width=0.2,alpha = 0.3)+
  scale_fill_manual(values = c("#f79f1f","#a3cb38","#1289a7"))+
  stat_compare_means(
    comparisons = compare_list1,
    method = "wilcox.test",
    label = "p.signif",
    size=6)+
  theme_bw()+
  theme(legend.position="none", panel.grid = element_blank(),
        axis.text.x=element_text(colour="black",size=16), #????x???潭缺?签??????????
        axis.text.y=element_text(size=16), #????y???潭缺?签??????????
        axis.title.y=element_text(size = 18), #????y???谋?????????????
        axis.title.x=element_text(size = 18),)+ 
  scale_y_continuous(expand = c(0.1, 0.1))+
  ylab ("Survival time without transfusion(months)")+
  xlab("(chr16:282411_G/A)\n(rs199594877)")

p1
ggsave(plot = p1, filename = "C://new_H_disease_NatureGenetics//Fig4-postGWAS//genotype-phenotype//case-control//HCS-chr16_282411_G_A_Survival-time-without-transfusion_pheno_snp_assoc.pdf", width = 6, height = 6)

