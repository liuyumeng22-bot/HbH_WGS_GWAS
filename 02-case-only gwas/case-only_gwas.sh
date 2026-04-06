source ~/.bashrc
gcta=/public/home/liuyumeng/software/GCTA/gcta-1.94.1-linux-kernel-3-x86_64/gcta64
path=/public/home/liuyumeng/new_gwas_results/585caseonly/plinkqc
outdir=/public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results
##计算GRM
$gcta --bfile $path/585.allchr_mind0.05_geno0.05_maf0.05_hwe --make-grm  --thread-num 10 --out $outdir/585case
########计算PCA
$gcta --grm $outdir/585case  --pca 10 --out $outdir/585case_pca


########计算稀松grm
$gcta --bfile $path/585.allchr_mind0.05_geno0.05_maf0.05_hwe --make-grm --sparse-cutoff 0.05 --thread-num 10 --out $outdir/585case_sp


########################GWAS:
source ~/.bashrc
gcta=/public/home/liuyumeng/software/GCTA/gcta-1.94.1-linux-kernel-3-x86_64/gcta64

bfilepath=/public/home/liuyumeng/new_gwas_results/585caseonly/plinkqc
path=/public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results
outdir=/public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results

for i in {1..10}
do
$gcta --bfile $bfilepath/585.allchr_mind0.05_geno0.05_maf0.05_hwe --grm-sparse $path/585case_sp --fastGWA-mlm --pheno $path/case_pheno_continous.txt --mpheno $i --qcovar $path/case_qcovar.txt  --covar $path/case_covar.txt  --thread-num 10 --out $outdir/$i
done



#################
cd ./sig.results

查看存在显著位点的GWAS结果文件：
for i in {1..51}; do awk '{if($10<8.54e-9) print $0}' ../$i.fastGWA >./$i.assoc.sig.txt; done
#####################
#################################
#################找独立independent SNP：（p<8.54e-9)
1          2           3          4         5        6        7          8           9       10
CHR     SNP     POS     A1      A2      N       AF1     BETA    SE      P

.ma文件：
SNP A1 A2 freq b se p N 

for i in {5,7,24,39,40,41,44,46,47,48}
 do awk '{print$2"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9"\t"$10"\t"$6}' /public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/$i.fastGWA >./$i.fast.ma
 done




######################是否是独立关联：
#####################cojo:
gcta=/public/home/liuyumeng/software/GCTA/gcta-1.94.1-linux-kernel-3-x86_64/gcta64

path1=/public/home/liuyumeng/new_gwas_results/585caseonly/plinkqc
outdir1=/public/home/liuyumeng/new_gwas_results/585caseonly/fast_mlm_results/assoc.results/sig.results/independent.loci

$gcta --bfile  $path1/585.allchr_mind0.05_geno0.05_maf0.05_hwe --cojo-file $outdir1/5.fast.ma --cojo-slct --cojo-p 5e-8 --out $outdir1/5pheno_5e-8_case

$gcta --bfile  $path1/585.allchr_mind0.05_geno0.05_maf0.05_hwe --cojo-file $outdir1/5.fast.ma --cojo-slct --cojo-p 8.57e-9 --out $outdir1/5pheno_8.57e-9_case



###################################################
#####################################
#############################
#####################GWAS catalog:
与LD信息结合查询：
Catalog本身不直接提供LD信息。您需要结合 LDlink 等工具，先确定您先导SNP的LD区域（例如，找出所有与您的SNP在目标人群中 r² > 0.8 的SNP），然后用这些SNP的rsID或区域去批量查询Catalog。


############open targets genetics:
1.known association（已经报道）
2.previously unidentified trait association(已知SNP,但是是新表型）
3.previously unidentified （全新为报道过的SNP,包括LD区域）