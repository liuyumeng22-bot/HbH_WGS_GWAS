#(文章材料方法书写流程：variants QC --Sample QC --VEP --Variants统计 --GWAS QC(PLINK QC))
#①VQSR(SNP:99.5;INDELS:99.0)
#②HF(DP>=10,GQ>=20),missing 10%,LongINDEL (length>20bp) were excluded(具体看小艾截图ppt)
#③Sample QC:sex check ....
#④VEP注释




#1.QC:GQ<20,DP<10,INDEL>20bp,missing>10%

source ~/.bashrc
conda activate bcftools
bcftools=/public/home/liuyumeng/anaconda3/envs/bcftools/bin/bcftools
path=/public/share/wchirdzhq2022/LiuYuMeng/591case_1246control.concatvcf
outdir=/public/share/wchirdzhq2022/LiuYuMeng/591case_1246control.concatvcf/QC/HF
######To assure genotypes with high confidence and accuracy, we further set individual genotypes with GQ<20 or DP<10 to missing data.
$bcftools filter $path/chr22.concat.vcf.gz -i 'FMT/GQ>=20 || FMT/DP>=10' --threads 10 -Oz -o $outdir/chr22.concat.GQDP.vcf.gz
$bcftools index -t $outdir/chr22.concat.GQDP.vcf.gz

######LongINDEL (length>20bp) were excluded
$bcftools filter --threads 5 -e 'ILEN>=20 || ILEN<=-20'  $outdir/chr22.concat.GQDP.vcf.gz -Oz -o $outdir/chr22.concat.GQDP.INDEL20.vcf.gz
$bcftools index -t $outdir/chr22.concat.GQDP.INDEL20.vcf.gz

###After such filtering, variants with more than 10% of missing genotypes accross individuals were eliminated for downstream analyses.
$bcftools filter -e "F_MISSING > 0.1" $outdir/chr22.concat.GQDP.INDEL20.vcf.gz -Oz -o $outdir/chr22.concat.GQDP.INDEL20.missing10.vcf.gz
$bcftools index -t $outdir/chr22.concat.GQDP.INDEL20.missing10.vcf.gz



#2.将SNPID的.改为chr:pos_ref/alt:
source ~/.bashrc
conda activate bcftools
bcftools=/public/home/liuyumeng/anaconda3/envs/bcftools/bin/bcftools

path=/public/share/wchirdzhq2022/LiuYuMeng/591case_1246control.concatvcf/QC/HF
outdir=/public/share/wchirdzhq2022/LiuYuMeng/591case_1246control.concatvcf/QC/change_snpid_vcf

######将SNPID的.改为chr:pos_ref/alt:
zcat $path/chr22.concat.GQDP.INDEL20.missing10.vcf.gz | awk -v OFS="\t" '/^#/|| $3=$1":"$2"_"$4"/"$5{print $0}' >$outdir/chr22.concat.GQDP.INDEL20.missing10.snpid.vcf

$bcftools view $outdir/chr22.concat.GQDP.INDEL20.missing10.snpid.vcf -Oz -o $outdir/chr22.concat.GQDP.INDEL20.missing10.snpid.vcf.gz

$bcftools index -t $outdir/chr22.concat.GQDP.INDEL20.missing10.snpid.vcf.gz



#3.###############剔除质控不合格样本
source ~/.bashrc
conda activate bcftools
bcftools=/public/home/liuyumeng/anaconda3/envs/bcftools/bin/bcftools

path=/public/share/wchirdzhq2022/LiuYuMeng/591case_1246control.concatvcf/QC/change_snpid_vcf
outdir1=/public/share/wchirdzhq2022/LiuYuMeng/591case_1246control.concatvcf/QC/sample_vcf/585HbHcase_vcf
outdir2=/public/share/wchirdzhq2022/LiuYuMeng/591case_1246control.concatvcf/QC/sample_vcf/585case_1246control_vcf

###########提取585H病样本vcf
$bcftools view -S /public/home/liuyumeng/585case.sample.list $path/chr22.concat.GQDP.INDEL20.missing10.snpid.vcf.gz -Oz -o $outdir1/chr22_585case.concat.GQDP.INDEL20.missing10.snpid.vcf.gz
$bcftools index -t $outdir1/chr22_585case.concat.GQDP.INDEL20.missing10.snpid.vcf.gz
$bcftools query -l $outdir1/chr22_585case.concat.GQDP.INDEL20.missing10.snpid.vcf.gz >/public/home/liuyumeng/sample.list/585case.sample.list

###########提取不合格样本：
$bcftools view -S ^/public/home/liuyumeng/remove.sample.list $path/chr22.concat.GQDP.INDEL20.missing10.snpid.vcf.gz -Oz -o $outdir2/chr22_585_1246.concat.GQDP.INDEL20.missing10.snpid.vcf.gz
$bcftools index -t $outdir2/chr22_585_1246.concat.GQDP.INDEL20.missing10.snpid.vcf.gz
$bcftools query -l $outdir2/chr22_585_1246.concat.GQDP.INDEL20.missing10.snpid.vcf.gz >/public/home/liuyumeng/sample.list/585_1246.sample.list


#4.###################585sample22个chr.vcf合并
source ~/.bashrc
conda activate bcftools
bcftools=/public/home/liuyumeng/anaconda3/envs/bcftools/bin/bcftools
path=/public/share/wchirdzhq2022/LiuYuMeng/591case_1246control.concatvcf/QC/sample_vcf/585HbHcase_vcf

ls $path/*.vcf.gz >$path/22chr.585.vcf.list
$bcftools concat -f $path/22chr.585.vcf.list -Oz -o $path/585.concat.allchr.vcf.gz
$bcftools index -t $path/585.concat.allchr.vcf.gz
