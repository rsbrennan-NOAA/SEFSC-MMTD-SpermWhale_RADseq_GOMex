### currentNe

#use only the major chromosomes, so I know the actual distance between snps.


module load bio/plink/1.90b6.23 

cd ~/spermWhaleRad/analysis/Ne

# need to add snp ids:

module load bio/bcftools/1.11
module load lib64/htslib/1.11

tabix ~/spermWhaleRad/analysis/Ne/freebayes_all_chr.vcf 
bcftools annotate ~/spermWhaleRad/analysis/Ne/freebayes_all_chr.vcf.gz -x ID -I +'%CHROM:%POS' > ~/spermWhaleRad/analysis/Ne/freebayes_all_chr_ids.vcf 

# run current Ne
~/spermWhaleRad/analysis/Ne/currentNe/currentNe ~/spermWhaleRad/analysis/Ne/freebayes_all_chr_ids.vcf 4 -o vcf_input
