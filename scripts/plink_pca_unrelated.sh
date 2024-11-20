#!/bin/bash
#SBATCH --job-name=plink_pca
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/spermWhaleRad/logout
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --partition=standard
#SBATCH --time=30:00

module load bio/plink/1.90b6.23

cd ~/spermWhaleRad/analysis/pop_structure/

#ld thin first, remove indivs, pca

plink --vcf ~/spermWhaleRad/analysis/freebayes/freebayes_ldthin_unrelated.vcf.gz \
--make-bed --out variants_unrelated_NoLD \
--allow-extra-chr --double-id

plink --bfile variants_unrelated_NoLD \
  --pca --out variants_NoLD_unrelated_PCA --allow-extra-chr --double-id

# remove indivs, ldthin, then pca
# full vcf

plink --vcf ~/spermWhaleRad/analysis/freebayes/freebayes_filtered_unrelated.vcf.gz \
--indep-pairwise 50 5 0.2 --allow-extra-chr --double-id \
--out variants_pruned_unrelated

# make thinned vcf
module load bio/vcftools/0.1.16

vcftools --gzvcf ~/spermWhaleRad/analysis/freebayes/freebayes_filtered_unrelated.vcf.gz --snps variants_pruned_unrelated.prune.in  --recode --recode-INFO-all --stdout | \ 
    bgzip > ~/spermWhaleRad/analysis/freebayes/freebayes_filtered_unrelated_ldthin.vcf.gz


plink --vcf ~/spermWhaleRad/analysis/freebayes/freebayes_filtered_unrelated_ldthin.vcf.gz \
--make-bed --out variants_unrelated_NoLD_2 \
--allow-extra-chr --double-id

plink --bfile variants_unrelated_NoLD_2 \
  --pca --out variants_unrelated_NoLD_2_PCA --allow-extra-chr --double-id
