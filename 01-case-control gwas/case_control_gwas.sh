#####################make grm:
source ~/.bashrc
gcta=/public/home/liuyumeng/software/GCTA/gcta-1.94.1-linux-kernel-3-x86_64/gcta64

path1=/public/home/liuyumeng/new_gwas_results/case-control/plinkQC/bfile
outdir1=/public/home/liuyumeng/new_gwas_results/case-control/result
 
  
 ###############
for i in {1..22}
do
echo $path1"/QC_chr"$i"_585_1246_geno0.05_maf0.05_hwe" >> $path1/bfile.chr.txt
done


########计算稀松grm
$gcta --mbfile $path1/bfile.chr.txt --make-grm --sparse-cutoff 0.05 --thread-num 10 --out $outdir1/new_case_control_sp
##计算GRM
$gcta --mbfile $path1/bfile.chr.txt --make-grm  --thread-num 10 --out $outdir1/new_case_control
########计算PCA
$gcta --grm $outdir1/new_case_control  --pca 10 --out $outdir1/new_pca


########################GWAS:
gcta=/public/home/liuyumeng/software/GCTA/gcta-1.94.1-linux-kernel-3-x86_64/gcta64

path=/public/home/liuyumeng/new_gwas_results/case-control/plinkQC/bfile
outdir=/public/home/liuyumeng/new_gwas_results/case-control/result

$gcta --mbfile $path/bfile.chr.txt --grm-sparse $outdir/new_case_control_sp --fastGWA-mlm-binary --pheno $outdir/585_1246sample.pheno.txt --qcovar $outdir/n    ew_case-control_qcovar.txt --covar $outdir/new_case-control_bcovar.txt   --thread-num 10 --out $outdir/new_case_control
~



#####################cojo:

awk '{print$2"\t"$4"\t"$5"\t"$7"\t"$11"\t"$12"\t"$13"\t"$6}' ./new_case_control.fastGWA >./conditional_analysis_independent/new_case_control.ma
#####################cojo:(8.57e-9)
##1.筛选出5e-8or 1e-5  的independent SNPs:
gcta=/public/home/liuyumeng/software/GCTA/gcta-1.94.1-linux-kernel-3-x86_64/gcta64

path1=/public/home/liuyumeng/new_gwas_results/case-control/plinkQC/bfile
outdir1=/public/home/liuyumeng/new_gwas_results/case-control/result/conditional_analysis_independent

$gcta --mbfile  $path1/bfile.chr.txt --cojo-file $outdir1/case_control.ma --cojo-slct --cojo-p 5e-8 --out $outdir1/new_case_control.ma



