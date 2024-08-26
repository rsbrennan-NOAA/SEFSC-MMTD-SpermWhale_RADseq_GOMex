#!/bin/bash
#SBATCH --job-name=filter2
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/spermWhaleRad/logout
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --time=30:00

source ~/.bashrc

mamba activate vcflib-1.0.9

module load bio/vcftools/0.1.16


vcftools --gzvcf ~/spermWhaleRad/analysis/freebayes/filtered.4.vcf.gz --remove-indv Pmac075 --recode --recode-INFO-all --stdout |  bgzip > ~/spermWhaleRad/analysis/freebayes/filtered.5.vcf.gz 

vcftools --gzvcf ~/spermWhaleRad/analysis/freebayes/filtered.5.vcf.gz  --max-missing 0.7 --maf 0.05 --min-meanDP 10 --recode --recode-INFO-all  --stdout >  ~/spermWhaleRad/analysis/freebayes/filtered.6.vcf

vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" ~/spermWhaleRad/analysis/freebayes/filtered.6.vcf |bgzip > ~/spermWhaleRad/analysis/freebayes/filtered.7.vcf.gz

#depths
vcftools --gzvcf ~/spermWhaleRad/analysis/freebayes/filtered.7.vcf.gz --site-mean-depth --out ~/spermWhaleRad/analysis/freebayes/filtered.7.depth


