################plink转换vcf文件为bfile格式
###############(1)vcf格式转换成plink格式
source ~/.bashrc
conda activate bcftools
bcftools=/public/home/liuyumeng/anaconda3/envs/bcftools/bin/bcftools
plink=/public/home/liuyumeng/software/plink1.9/plink

path=/public/home/liuyumeng/RVS_VCF
outdir=/public/home/liuyumeng/new_gwas_results/585caseonly/plinkqc

$plink --vcf $path/maf_0_filter_585concat.GQDP.INDEL20.missing10.vcf.gz --allow-extra-chr  --const-fid 1 --threads 10  --out $outdir/585.allchr

####质控：
$plink --bfile $outdir/585.allchr  --mind 0.05 --make-bed --threads 10  --out $outdir/585.allchr_mind0.05

$plink --bfile $outdir/585.allchr_mind0.05 --geno 0.05 --maf 0.05  --hwe 1e-6 --make-bed --threads 10  --out $outdir/585.allchr_mind0.05_geno0.05_maf0.05_hwe
