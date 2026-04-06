###################提取相应列画manhattan,qq图：
cd /public/home/liuyumeng/new_gwas_results/case-control/result/manhattan.plot
awk '{print$2"\t"$1"\t"$3"\t"$13}' ../new_case_control.fastGWA>./assoc.manhatant.qq.txt
awk '!/chr21:5219327_TA\/T/' ./assoc.manhatant.qq.txt >./delete_1_assoc.manhatant.qq.txt

##################提取相应列和P<1e-5独立位点用于后续VEP注释
source ~/.bashrc
conda activate bcftools
bcftools=/public/home/liuyumeng/anaconda3/envs/bcftools/bin/bcftools
path1=/public/home/liuyumeng/new_gwas_results/case-control/result/vep/chr_pos
path2=/public/home/liuyumeng/591case_1246control.concatvcf/QC/sample_vcf/585case_1246control_vcf
outdir=/public/home/liuyumeng/new_gwas_results/case-control/result/vep
$bcftools view -R $path1/chr1.txt $path2/chr1_585_1246.concat.GQDP.INDEL20.missing10.snpid.vcf.gz -Oz -o $path1/chr1.vcf.gz
$bcftools view -R $path1/chr2.txt $path2/chr2_585_1246.concat.GQDP.INDEL20.missing10.snpid.vcf.gz -Oz -o $path1/chr2.vcf.gz
$bcftools view -R $path1/chr5.txt $path2/chr5_585_1246.concat.GQDP.INDEL20.missing10.snpid.vcf.gz -Oz -o $path1/chr5.vcf.gz
$bcftools view -R $path1/chr6.txt $path2/chr6_585_1246.concat.GQDP.INDEL20.missing10.snpid.vcf.gz -Oz -o $path1/chr6.vcf.gz
$bcftools view -R $path1/chr9.txt $path2/chr9_585_1246.concat.GQDP.INDEL20.missing10.snpid.vcf.gz -Oz -o $path1/chr9.vcf.gz
$bcftools view -R $path1/chr12.txt $path2/chr12_585_1246.concat.GQDP.INDEL20.missing10.snpid.vcf.gz -Oz -o $path1/chr12.vcf.gz
$bcftools view -R $path1/chr16.txt $path2/chr16_585_1246.concat.GQDP.INDEL20.missing10.snpid.vcf.gz -Oz -o $path1/chr16.vcf.gz
$bcftools view -R $path1/chr17.txt $path2/chr17_585_1246.concat.GQDP.INDEL20.missing10.snpid.vcf.gz -Oz -o $path1/chr17.vcf.gz
$bcftools view -R $path1/chr20.txt $path2/chr20_585_1246.concat.GQDP.INDEL20.missing10.snpid.vcf.gz -Oz -o $path1/chr20.vcf.gz
$bcftools view -R $path1/chr21.txt $path2/chr21_585_1246.concat.GQDP.INDEL20.missing10.snpid.vcf.gz -Oz -o $path1/chr21.vcf.gz
for i in $path1/*.vcf.gz
do
$bcftools index -t $i
done
$bcftools concat $path1/chr1.vcf.gz $path1/chr2.vcf.gz $path1/chr5.vcf.gz $path1/chr6.vcf.gz  $path1/chr9.vcf.gz  $path1/chr12.vcf.gz  $path1/chr16.vcf.gz $path1/chr17.vcf.gz $path1/chr20.vcf.gz $path1/chr21.vcf.gz -Oz -o $outdir/chr.vep.vcf.gz
$bcftools index -t $outdir/chr.vep.vcf.gz

#############
###################################VEP注释：
source ~/.bashrc
conda activate vep
reference=/public/home/liuyumeng/reference
path=/public/home/liuyumeng/new_gwas_results/case-control/result/vep
vep -i $path/chr.vep.vcf.gz  --assembly GRCh38 --cache --cache_version 108 --dir /public/home/liuyumeng/vep/reference/  --port 3337 --offline --fasta $reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna  --everything -o $path/vep.txt

########################
####################################locuszooms:
#三个leading SNPS:
chr16:336749_G/A
chr20:55443626_AGGAGTGGGGCCT/A
chr2:225873608_TTG/T
chr16:282411_G/A
path=/public/home/liuyumeng/new_gwas_results/case-control/plinkQC/bfile
outdir=/public/home/liuyumeng/new_gwas_results/case-control/result/locuszoom
plink=/public/home/liuyumeng/software/plink1.9/plink
# 使用PLINK生成LD文件
$plink --bfile $path/QC_chr16_585_1246_geno0.05_maf0.05_hwe --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --r2 --out $outdir/chr16_ld_results
$plink --bfile $path/QC_chr2_585_1246_geno0.05_maf0.05_hwe --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --r2 --out $outdir/chr2_ld_results
$plink --bfile $path/QC_chr20_585_1246_geno0.05_maf0.05_hwe --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --r2 --out $outdir/chr20_ld_results
 # 提取并格式化LD文件
awk 'BEGIN {print "snp1 snp2 dprime rsquare"} {print $3, $6, "NA", $7}' $outdir/chr16_ld_results.ld > $outdir/chr16_ld_file.txt
awk 'BEGIN {print "snp1 snp2 dprime rsquare"} {print $3, $6, "NA", $7}' $outdir/chr2_ld_results.ld > $outdir/chr2_ld_file.txt
awk 'BEGIN {print "snp1 snp2 dprime rsquare"} {print $3, $6, "NA", $7}' $outdir/chr20_ld_results.ld > $outdir/chr20_ld_file.txt

cd /public/home/liuyumeng/new_gwas_results/case-control/result/manhattan.plot
awk 'BEGIN {print "MarkerName\tP-value"} {print "chr"$2":"$3"\t"$4}' ./assoc.manhatant.qq.txt >../locuszoom/case_control.metal
###locuszoom:
source ~/.bashrc
conda activate locuszoom
plink=/public/home/liuyumeng/software/plink1.9/plink
locuszoom=/public/home/liuyumeng/software/locuszoom/bin/locuszoom
path=/public/home/liuyumeng/new_gwas_results/case-control/result/locuszoom
outdir=/public/home/liuyumeng/new_gwas_results/case-control/result/locuszoom
python2 $locuszoom --metal $path/case_control.metal --refsnp chr16:336749_G/A --chr 16 --start 34870 --end 1034870 --build hg38 --ld $outdir/chr16_ld_file.txt --ignore-vcf-filter  --prefix $outdir/chr16_336749_G_A    --plotonly
chr16:282411_G/A
python2 $locuszoom --metal $path/case_control.metal --refsnp chr16:282411_G/A --chr 16 --start 34870 --end 1034870 --build hg38 --ld $outdir/chr16_ld_file.txt --ignore-vcf-filter  --prefix $outdir/chr16_282411_G_A    --plotonly
python2 $locuszoom --metal $path/case_control.metal --refsnp chr20:55443626_AGGAGTGGGGCCT/A --chr 20 --start 54943626 --end 55943626 --build hg38 --ld $outdir/chr20_ld_file.txt --ignore-vcf-filter  --prefix $outdir/chr20_55443626_AGGAGTGGGGCCT_A    --plotonly
python2 $locuszoom --metal $path/case_control.metal --refsnp chr2:225873608_TTG/T --chr 2 --start 225373608 --end 226373608 --build hg38 --ld $outdir/chr2_ld_file.txt --ignore-vcf-filter  --prefix $outdir/chr2_225873608_TTG_T   --plotonly




###########
###############################finemap/SUSIER
#寻找lead snp: chr16:282411_G/A
#以lead snp为中心，上下游500kb的范围的所有snp列表，制成z.file
#chr16:282411_G/A
cd /public/home/liuyumeng/new_gwas_results/case-control/result/SusieR
awk -v c="16" -v p="282411" -v w="500000" 'NR==1 || ($1==c && $3>=p-w && $3<=p+w)' ../new_case_control.fastGWA >./chr16_282411_G_A.txt
awk  'NR==1 || ($1==16 && $3>=34870 && $3<=1034870)' ../new_case_control.fastGWA >./chr16_282411_G_A.txt

# 跳过标题行，计算zscore = BETA(第11列)/SE(第12列)
awk 'NR==1 {next} {zscore=$11/$12; print $2, zscore, $6}' OFS='\t' ./chr16_282411_G_A.txt >./chr16_282411_G_A.zscore.txt

#制作SUSIER
source ~/.bashrc
conda activate bcftools
bcftools=/public/home/liuyumeng/anaconda3/envs/bcftools/bin/bcftools
plink=/public/home/liuyumeng/software/plink1.9/plink

path=/public/home/liuyumeng/new_gwas_results/case-control/plinkQC/bfile
outdir=/public/home/liuyumeng/new_gwas_results/case-control/result/SusieR

$plink --bfile  $path/QC_chr16_585_1246_geno0.05_maf0.05_hwe --extract $outdir/chr16.snp.list --make-bed --out $outdir/chr16

$plink --bfile $outdir/chr16 --r2 --matrix --out $outdir/chr16.LD.matrix


> summary(fitted_rss3)$cs
   cs cs_log10bf cs_avg_r2 cs_min_r2 variable
1   1  135.06033 1.0000000 1.0000000      160
2   2  231.13085 1.0000000 1.0000000      217
3   3  252.76173 1.0000000 1.0000000      214
4   4   94.33804 1.0000000 1.0000000      276
5   5   74.47540 1.0000000 1.0000000      102
6   6   88.20157 1.0000000 1.0000000      748
7   8   94.63374 1.0000000 1.0000000      516
8   9   70.12614 1.0000000 1.0000000      564
9  10   60.96306 1.0000000 1.0000000      136
10  7   78.02365 0.9623021 0.9623021  265,268
> top_snps <- c(160, 217, 214, 276, 102, 748, 516, 564, 136, 265, 268)
> # 确认这些SNP的PIP是否接近1
> top_pips <- fitted_rss3$pip[top_snps]
> data.frame(SNP=top_snps, PIP=top_pips)
   SNP        PIP
1  160 1.00000000
2  217 1.00000000
3  214 1.00000000
4  276 0.96766208
5  102 1.00000000
6  748 0.99999999
7  516 1.00000000
8  564 1.00000000
9  136 0.99999818
10 265 0.07280176
11 268 0.92600415
> eq[top_snps, ]  
                   SNP   zscore    N
160   chr16:282411_G/A  20.1655 1831
217   chr16:342055_T/A -19.4931 1741
214   chr16:336749_G/A  18.6835 1820
276   chr16:391710_T/C -17.4277 1830
102   chr16:159817_T/C -19.2805 1814
748   chr16:663726_C/G  18.3713 1831
516 chr16:517249_TCA/T  19.9277 1831
564   chr16:555481_T/C -15.0242 1830
136   chr16:247751_T/C -16.9848 1831
265   chr16:389618_G/T -16.9923 1831
268   chr16:389821_G/T -17.1823 1829



#################################LDBlockShow图： 
awk '{print"chr"$1"\t"$3"\t"$13}' ./case_control.fastGWA >./case_control.gwas.pvalue


LDBlockShow=/public/home/liuyumeng/software/LDBlockShow-1.40/bin/LDBlockShow
path1=/public/share/wchirdzhq2022/LiuYuMeng/591case_1246control.concatvcf/QC/sample_vcf/585case_1246control_vcf
path2=/public/home/liuyumeng/GWAS.results/case_control/finalresult/LD_prune
outdir=/public/home/liuyumeng/GWAS.results/case_control/finalresult/LD_prune/LDBlockShow

$LDBlockShow -i $path1/chr16_change_sampleid_585_1246.concat.GQDP.INDEL20.missing10.snpid.vcf.gz  -o $outdir/chr16_342055 -InGWAS $path2/case_control.gwas.pvalue -r chr16:34870:1032254 -OutPng -SeleVar 3 -TopSite


###################################################genotype-phenotype 
chr16:266778_CAT/C
chr16:282411_G/A
chr16:337037_A/C
chr16:342055_T/A

plink=/public/home/liuyumeng/software/plink1.9/plink

path=/public/share/wchirdzhq2022/LiuYuMeng/591case_1246control.concatvcf/QC/bfile/case_control/bfile/LDprune_0.9
outdir=/public/home/liuyumeng/GWAS.results/case_control/finalresult/new_LD_prune/result/genotyep_phenotype

$plink --bfile $path/chr16_LDprune_0.9_pruned --extract $outdir/chr16.causal.snp.list --recode --out $outdir/chr16.causal.snp
