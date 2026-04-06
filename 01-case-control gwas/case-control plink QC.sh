################plink转换vcf文件为bfile格式
############case-control:
#########vcf格式转换成plink格式
source ~/.bashrc
conda activate bcftools
bcftools=/public/home/liuyumeng/anaconda3/envs/bcftools/bin/bcftools
plink=/public/home/liuyumeng/software/plink1.9/plink


outdir1=/public/home/liuyumeng/591case_1246control.concatvcf/QC/bfile/case_control/bfile
outdir2=/public/home/liuyumeng/new_gwas_results/case-control/plinkQC/bfile
outdir3=/public/home/liuyumeng/new_gwas_results/case-control/plinkQC/bfile/geno_0.05_maf_0.05_hwe_1e-6_LDprune0.9

#########质控：

$plink --bfile $outdir1/chr22_585_1246 --geno 0.05 --maf 0.05  --hwe 1e-6 --make-bed --out $outdir2/QC_chr22_585_1246_geno0.05_maf0.05_hwe --threads 10

$plink --bfile $outdir2/QC_chr22_585_1246_geno0.05_maf0.05_hwe --indep-pairwise 1000 100 0.9  --out $outdir3/QC_chr22_585_1246_geno0.05_maf0.05_hwe_LD_prune0.9 
--threads 10

$plink --bfile $outdir2/QC_chr22_585_1246_geno0.05_maf0.05_hwe  --extract  $outdir3/QC_chr22_585_1246_geno0.05_maf0.05_hwe_LD_prune0.9.prune.in  --make-bed --ou
t $outdir3/QC_chr22_585_1246_geno0.05_maf0.05_hwe_LD_pruned  --threads 10




