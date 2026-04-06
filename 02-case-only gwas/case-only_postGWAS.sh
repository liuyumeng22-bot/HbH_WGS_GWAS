##########################
查看存在显著位点的GWAS结果文件：
cd /public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/sig.results
for i in {1..51}; do awk '{if($10<8.54e-9) print $0}' ../$i.fastGWA >./$i.assoc.sig.txt; done
提取想要列用于后续画qq plot,manhattan plot：
awk '{print$2"\t"$1"\t"$3"\t"$10"\t""TBIL"}' ../5.fastGWA >./5.assoc.txt
 awk '{print$2"\t"$1"\t"$3"\t"$10"\t""IBIL"}' ../7.fastGWA >./7.assoc.txt
 awk '{print$2"\t"$1"\t"$3"\t"$10"\t""BAS2"}' ../24.fastGWA >./24.assoc.txt
 awk '{print$2"\t"$1"\t"$3"\t"$10"\t""HbA2"}' ../39.fastGWA >./39.assoc.txt
 awk '{print$2"\t"$1"\t"$3"\t"$10"\t""HbCS2"}' ../40.fastGWA >./40.assoc.txt
 awk '{print$2"\t"$1"\t"$3"\t"$10"\t""HbCS"}' ../41.fastGWA >./41.assoc.txt
 awk '{print$2"\t"$1"\t"$3"\t"$10"\t""HbBarts"}' ../44.fastGWA >./44.assoc.txt
 awk '{print$2"\t"$1"\t"$3"\t"$10"\t""TIBC"}' ../46.fastGWA >./46.assoc.txt
 awk '{print$2"\t"$1"\t"$3"\t"$10"\t""UIBC"}' ../47.fastGWA >./47.assoc.txt
 awk '{print$2"\t"$1"\t"$3"\t"$10"\t""TFR"}' ../48.fastGWA >./48.assoc.txt
将各个表型的文件合并
cat ./*.assoc.txt >merge.assoc.txt
awk '{if($4<=0.05) print $0}' ./merge.assoc.txt >./merge.assoc.p0.05.txt 
cat ./*.sig.txt >./merge.sig.txt
awk '{print"chr"$1"\t"$3}' ./merge.sig.txt >./merge.sig.snp.txt

###############################提取leadsnp的vcf用于以后的vep
source ~/.bashrc
conda activate bcftools
bcftools=/public/home/liuyumeng/anaconda3/envs/bcftools/bin/bcftools
path1=/public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/sig.results/independent.loci/vep
path2=/public/home/liuyumeng/RVS_VCF
outdir=/public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/sig.results/independent.loci/vep
$bcftools view -R $path1/merged.SNP.POS.txt $path2/maf_0_filter_585concat.GQDP.INDEL20.missing10.vcf.gz -Oz -o $outdir/vep.vcf.gz
$bcftools index -t $outdir/vep.vcf.gz

###################################VEP注释：
source ~/.bashrc
conda activate vep
reference=/public/home/liuyumeng/reference
path=/public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/sig.results/independent.loci/vep
vep -i $path/vep.vcf.gz  --assembly GRCh38 --cache --cache_version 108 --dir /public/home/liuyumeng/vep/reference/  --port 3337 --offline --fasta $reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna  --everything -o $path/vep.txt


#############################locuszoom
awk 'BEGIN {print "MarkerName\tP-value"} {print "chr"$2":"$3"\t"$4}' ../5.assoc.txt > ./TBIL.metal
awk 'BEGIN {print "MarkerName\tP-value"} {print "chr"$2":"$3"\t"$4}' ../7.assoc.txt > ./IBIL.metal
awk 'BEGIN {print "MarkerName\tP-value"} {print "chr"$2":"$3"\t"$4}' ../24.assoc.txt >./BAS2.metal
awk 'BEGIN {print "MarkerName\tP-value"} {print "chr"$2":"$3"\t"$4}'  ../39.assoc.txt >./HbA2.metal
awk 'BEGIN {print "MarkerName\tP-value"} {print "chr"$2":"$3"\t"$4}'  ../41.assoc.txt >./HbCS.metal
awk 'BEGIN {print "MarkerName\tP-value"} {print "chr"$2":"$3"\t"$4}' ../44.assoc.txt >./HbBarts.metal
awk 'BEGIN {print "MarkerName\tP-value"} {print "chr"$2":"$3"\t"$4}'  ../46.assoc.txt >./TIBC.metal
awk 'BEGIN {print "MarkerName\tP-value"} {print "chr"$2":"$3"\t"$4}'  ../47.assoc.txt >./UIBC.metal
awk 'BEGIN {print "MarkerName\tP-value"} {print "chr"$2":"$3"\t"$4}' ../48.assoc.txt >./TRF.metal
source ~/.bashrc
conda activate locuszoom
plink=/public/home/liuyumeng/software/plink1.9/plink
locuszoom=/public/home/liuyumeng/software/locuszoom/bin/locuszoom
path=/public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/sig.results/locuszoom
python2 $locuszoom --metal $path/TBIL.metal --refsnp chr2:233762816_G/T --flank 500kb --build hg38 --ld $path/case_only_ld_file.txt --ignore-vcf-filter  --prefix $path/chr2_233762816_G_T_TBIL    --plotonly
python2 $locuszoom --metal $path/BAS2.metal --refsnp chr2:231526555_T/C --flank 500kb --build hg38 --ld $path/case_only_ld_file.txt --ignore-vcf-filter --prefix $path/chr2_231526555_T_C_BAS2    --plotonly
python2 $locuszoom --metal $path/HbA2.metal --refsnp chr5:83721739_G/GCGCGCGCGCACACA --flank 500kb --build hg38 --ld $path/case_only_ld_file.txt --ignore-vcf-filter  --prefix $path/chr5_83721739_G_GCGCGCGCGCACACA_HbA2    --plotonly
python2 $locuszoom --metal $path/HbCS.metal --refsnp chr16:247672_C/A --chr 16 --start 34870 --end 1034870 --build hg38 --ld $path/case_only_ld_file.txt --ignore-vcf-filter  --prefix $path/chr16_247672_C_A_HbCS    --plotonly
python2 $locuszoom --metal $path/HbBarts.metal --refsnp chr6:135101158_A/T --flank 500kb --build hg38 --ld $path/case_only_ld_file.txt --ignore-vcf-filter  --prefix $path/chr6_135101158_A_T_HbBarts    --plotonly
python2 $locuszoom --metal $path/TIBC.metal --refsnp chr3:133760356_G/C --flank 500kb --build hg38 --ld $path/case_only_ld_file.txt --ignore-vcf-filter  --prefix $path/chr3_133760356_G_C_TIBC    --plotonly
python2 $locuszoom --metal $path/UIBC.metal --refsnp chr16:262254_A/G --chr 16 --start 34870 --end 1034870 --build hg38 --ld $path/case_only_ld_file.txt --ignore-vcf-filter  --prefix $path/chr16_262254_A_G_UIBC    --plotonly
python2 $locuszoom --metal $path/TRF.metal --refsnp chr3:133762874_C/T --flank 500kb --build hg38 --ld $path/case_only_ld_file.txt --ignore-vcf-filter  --prefix $path/chr3_133762874_C_T_TRF    --plotonly
#############################################
五. ###############################finemap/SUSIER
寻找lead snp:
以lead snp为中心，上下游500kb的范围的所有snp列表，制成z.file
Trait	SNP
BAS%	chr2:231526555_T/C √
HbA2	chr5:83721739_G/GCGCGCGCGCACACA
HbCS%	chr16:247672_C/A
HbCS	chr16:247672_C/A √
HbBarts	chr6:135101158_A/T √
TIBC	chr3:133760356_G/C
UIBC	chr16:262254_A/G √
TRF	chr3:133762874_C/T
TBIL	chr2:233762816_G/T
IBIL	chr2:233762816_G/T

SNP	zscore	N
BAS%	chr2:231526555_T/C sce qtl:NMUR1 CD8 Tcell
cd /public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/sig.results/SusieR
CHR     SNP     POS     A1      A2      N       AF1     BETA    SE      P
awk  'NR==1 || ($1==2 && $3>=231026555 && $3<=232026555)' /public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/24.fastGWA >./BAS2_chr2_231526555_T_C.txt
# 跳过标题行，计算zscore = BETA(第8列)/SE(第9列)
CHR     SNP     POS     A1      A2      N       AF1     BETA    SE      P
awk 'NR==1 {next} {zscore=$8/$9; print $2, zscore, $6}' OFS='\t' ./BAS2_chr2_231526555_T_C.txt >./BAS2_chr2_231526555_T_C.zscore.txt
################plink转换vcf文件为bfile格式
###############(1)vcf格式转换成plink格式
source ~/.bashrc
conda activate bcftools
bcftools=/public/home/liuyumeng/anaconda3/envs/bcftools/bin/bcftools
plink=/public/home/liuyumeng/software/plink1.9/plink
path=/public/home/liuyumeng/new_gwas_results/585caseonly/plinkqc
outdir=/public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/sig.results/SusieR
$plink --bfile  $path/585.allchr_mind0.05_geno0.05_maf0.05_hwe --extract $outdir/BAS2_chr2_231526555_T_C.SNP.txt --make-bed --out $outdir/BAS2_chr2_231526555_T_C --threads 10
$plink --bfile $outdir/BAS2_chr2_231526555_T_C --r2 --matrix --out $outdir/BAS2_chr2_231526555_T_C.LD.matrix --threads 10
  cs cs_log10bf cs_avg_r2 cs_min_r2 variable
1  1  102804.12         1         1      963
2  2   92422.48         1         1      966
3  3   18911.09         1         1      983
4  4   38201.26         1         1      976
> top_snps <- c(963, 966, 983, 976)
> # 确认这些SNP的PIP是否接近1
> top_pips <- fitted_rss3$pip[top_snps]
> data.frame(SNP=top_snps, PIP=top_pips)
  SNP PIP
1 963   1
2 966   1
3 983   1
4 976   1
> eq[top_snps, ]  
                   SNP   zscore   N
963 chr2:231519173_A/C -5.75641 570
966 chr2:231520664_C/G -5.83704 571
983 chr2:231523671_T/A  5.55632 571
976 chr2:231522280_C/T -5.47613 570


###########2.HbBarts	chr6:135101158_A/T √
SNP	zscore	N
BAS%	chr2:231526555_T/C sce qtl:NMUR1 CD8 Tcell
cd /public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/sig.results/SusieR/HbBarts_chr6_135101158_A_T
CHR     SNP     POS     A1      A2      N       AF1     BETA    SE      P
awk  'NR==1 || ($1==6 && $3>=134601158 && $3<=135601158)' /public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/44.fastGWA >./HbBarts_chr6_135101158_A_T.txt
# 跳过标题行，计算zscore = BETA(第8列)/SE(第9列)
CHR     SNP     POS     A1      A2      N       AF1     BETA    SE      P
awk 'NR==1 {next} {zscore=$8/$9; print $2, zscore, $6}' OFS='\t' ./HbBarts_chr6_135101158_A_T.txt >./HbBarts_chr6_135101158_A_T.zscore.txt
################plink转换vcf文件为bfile格式
###############(1)vcf格式转换成plink格式
source ~/.bashrc
conda activate bcftools
bcftools=/public/home/liuyumeng/anaconda3/envs/bcftools/bin/bcftools
plink=/public/home/liuyumeng/software/plink1.9/plink
path=/public/home/liuyumeng/new_gwas_results/585caseonly/plinkqc
outdir=/public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/sig.results/SusieR/HbBarts_chr6_135101158_A_T
$plink --bfile  $path/585.allchr_mind0.05_geno0.05_maf0.05_hwe --extract $outdir/HbBarts_chr6_135101158_A_T.snp.txt --make-bed --out $outdir/HbBarts_chr6_135101158_A_T --threads 10
$plink --bfile $outdir/HbBarts_chr6_135101158_A_T --r2 --matrix --out $outdir/HbBarts_chr6_135101158_A_T.LD.matrix --threads 10
##########SusieR
summary(fitted_rss3)$cs
           SNP         PIP
1    828 1.000000000
2   1087 1.000000000
3   1106 0.984507538
4   1058 0.999964297

                                  SNP   zscore   N
828                chr6:134962243_TA/T -3.48210 582
1087                chr6:135094070_G/A  5.52952 582
1106                chr6:135102274_G/A -2.45694 582
1058                chr6:135070962_C/A -2.01410 582



##################33UIBC_chr16_262254_A_G 
cd /public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/sig.results/SusieR/UIBC_chr16_262254_A_G
SNP	zscore	N
UIBC_chr16_262254_A_G
CHR     SNP     POS     A1      A2      N       AF1     BETA    SE      P
awk  'NR==1 || ($1==16 && $3>=34870&& $3<=1034870)' /public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/47.fastGWA >./UIBC_chr16_262254_A_G.txt
# 跳过标题行，计算zscore = BETA(第8列)/SE(第9列)
CHR     SNP     POS     A1      A2      N       AF1     BETA    SE      P
awk 'NR==1 {next} {zscore=$8/$9; print $2, zscore, $6}' OFS='\t' ./UIBC_chr16_262254_A_G.txt >./UIBC_chr16_262254_A_G.zscore.txt
################plink转换vcf文件为bfile格式
###############(1)vcf格式转换成plink格式
ssource ~/.bashrc
conda activate bcftools
bcftools=/public/home/liuyumeng/anaconda3/envs/bcftools/bin/bcftools
plink=/public/home/liuyumeng/software/plink1.9/plink
path=/public/home/liuyumeng/new_gwas_results/585caseonly/plinkqc
outdir=/public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/sig.results/SusieR/UIBC_chr16_262254_A_G
$plink --bfile  $path/585.allchr_mind0.05_geno0.05_maf0.05_hwe --extract $outdir/UIBC_chr16_262254_A_G.snp.txt --make-bed --out $outdir/UIBC_chr16_262254_A_G --threads 10
$plink --bfile $outdir/UIBC_chr16_262254_A_G --r2 --matrix --out $outdir/UIBC_chr16_262254_A_G.LD.matrix --threads 10
###########3summary(fitted_rss3)$cs
  cs cs_log10bf cs_avg_r2 cs_min_r2                                     variable
1  1   3.785116 0.8098285   0.60246 76,77,78,79,82,83,84,85,90,91,93,94,95,97,98
> top_snps <- c(76,77,78,79,82,83,84,85,90,91,93,94,95,97,98)
> # 确认这些SNP的PIP是否接近1
> top_pips <- fitted_rss3$pip[top_snps]
> data.frame(SNP=top_snps, PIP=top_pips)
   SNP        PIP
1   76 0.04841840
2   77 0.07213477
3   78 0.07213477
4   79 0.07213477
5   82 0.08261457
6   83 0.07643837
7   84 0.07498098
8   85 0.03835626
9   90 0.07469151
10  91 0.07469151
11  93 0.04581487
12  94 0.08846308
13  95 0.07154832
14  97 0.08825100
15  98 0.06890754
> eq[top_snps, ]  
                  SNP  zscore   N
76   chr16:244750_A/G 5.67714 567
77   chr16:245719_A/G 5.76510 579
78   chr16:245796_T/C 5.76510 579
79   chr16:246185_C/A 5.76510 579
82   chr16:250178_A/G 5.79402 579
83   chr16:250389_T/C 5.77756 579
84   chr16:250400_T/C 5.77406 579
85   chr16:250642_G/A 5.62636 579
90   chr16:255269_T/C 5.77379 579
91   chr16:255378_G/C 5.77379 579
93   chr16:261854_A/G 5.66537 579
94   chr16:262254_A/G 5.80858 579
95   chr16:263407_A/G 5.76385 579
97 chr16:268155_T/TCA 5.80731 579
98   chr16:268475_C/G 5.75474 579
                                                                               



##################4HbCS_chr16_247672_C_A 
cd /public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/sig.results/SusieR/HbCS_chr16_247672_C_A
SNP	zscore	N
HbCS_chr16_247672_C_A
CHR     SNP     POS     A1      A2      N       AF1     BETA    SE      P
awk  'NR==1 || ($1==16 && $3>=34870&& $3<=1034870)' /public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/41.fastGWA >./HbCS_chr16_247672_C_A.txt
# 跳过标题行，计算zscore = BETA(第8列)/SE(第9列)
CHR     SNP     POS     A1      A2      N       AF1     BETA    SE      P
awk 'NR==1 {next} {zscore=$8/$9; print $2, zscore, $6}' OFS='\t' ./HbCS_chr16_247672_C_A.txt >./HbCS_chr16_247672_C_A.zscore.txt
awk '{print$1}' ./HbCS_chr16_247672_C_A.zscore.txt >./HbCS_chr16_247672_C_A.snp.txt
################plink转换vcf文件为bfile格式
###############(1)vcf格式转换成plink格式
ssource ~/.bashrc
conda activate bcftools
bcftools=/public/home/liuyumeng/anaconda3/envs/bcftools/bin/bcftools
plink=/public/home/liuyumeng/software/plink1.9/plink
path=/public/home/liuyumeng/new_gwas_results/585caseonly/plinkqc
outdir=/public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/sig.results/SusieR/HbCS_chr16_247672_C_A
$plink --bfile  $path/585.allchr_mind0.05_geno0.05_maf0.05_hwe --extract $outdir/HbCS_chr16_247672_C_A.snp.txt--make-bed --out $outdir/HbCS_chr16_247672_C_A --threads 10
$plink --bfile $outdir/HbCS_chr16_247672_C_A --r2 --matrix --out $outdir/HbCS_chr16_247672_C_A.LD.matrix --threads 10
  summary(fitted_rss3)$cs
NULL                                                                                                                                                                                                                                                                                                                                                                       
################genotype_phenotype
plink=/public/home/liuyumeng/software/plink1.9/plink
path=/public/share/wchirdzhq2022/LiuYuMeng/591case_1246control.concatvcf/QC/bfile/585case/bfile/LD_prune0.9
outdir=/public/home/liuyumeng/GWAS.results/585case/final.result/new_LD_prune_GCTA/zscore/sig.result/genotype_phenotype
$plink --bfile $path/0.9LD_pruned --extract $outdir/all_leading.snp.list --recode --out $outdir/all_leading_snp
$plink --bfile $outdir1/score_chr16_241645_CT --r2 'chr16:254804_A/T' 'chr16:241645_C/T'
