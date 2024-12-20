#!/bin/bash
#SBATCH --job-name=filter3
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/spermWhaleRad/logout
#SBATCH -n 1
#SBATCH --mem=6G
#SBATCH --partition=standard
#SBATCH --time=30:00

source ~/.bashrc

mamba activate vcflib-1.0.9
module load bio/vcftools/0.1.16

# full vcf
vcftools --gzvcf ~/spermWhaleRad/analysis/freebayes/filtered.final_ids.vcf.gz --remove ~/spermWhaleRad/analysis/relatedness/relatedindivs.txt  --recode --recode-INFO-all --stdout  | bgzip > ~/spermWhaleRad/analysis/freebayes/freebayes_filtered_unrelated.vcf.gz

# ld thinned:
vcftools --gzvcf ~/spermWhaleRad/analysis/freebayes/filtered.final_LDthin.vcf.gz --remove ~/spermWhaleRad/analysis/relatedness/relatedindivs.txt  --recode --recode-INFO-all --stdout | bgzip   > ~/spermWhaleRad/analysis/freebayes/freebayes_ldthin_unrelated.vcf.gz

# 
