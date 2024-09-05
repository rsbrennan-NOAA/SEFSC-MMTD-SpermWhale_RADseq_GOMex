#!/bin/bash
#SBATCH --job-name=filter1
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/spermWhaleRad/logout
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --time=3:00:00

source ~/.bashrc

mamba activate vcflib-1.0.9

module load bio/vcftools
module load bio/bcftools

cd ~/spermWhaleRad/analysis/freebayes

# create a filtered VCF containing only invariant sites
#vcftools --gzvcf freebayes_allsites.vcf.gz \
#	--max-maf 0.05 \
#       --remove-indv Pmac075 --remove-indels \
#        --max-missing 0.7 \
#        --min-meanDP 10 \
#        --max-meanDP 80 \
#	--recode --stdout | bgzip -c > filtered_invariant.vcf.gz

# index both vcfs using tabix
tabix -f filtered_invariant.vcf.gz
tabix -f filtered.final.vcf.gz

# combine the two VCFs using bcftools concat
bcftools concat \
--allow-overlaps \
filtered.final.vcf.gz filtered_invariant.vcf.gz \
-O z -o combined_filtered_invariant.vcf.gz
