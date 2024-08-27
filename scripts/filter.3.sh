#!/bin/bash
#SBATCH --job-name=filter3
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

vcftools --gzvcf ~/spermWhaleRad/analysis/freebayes/filtered.7.vcf.gz --max-meanDP 53.4 --exclude-positions ~/spermWhaleRad/analysis/freebayes/HD_exclude.txt --recode --recode-INFO-all --stdout |  bgzip > ~/spermWhaleRad/analysis/freebayes/filtered.final.vcf.gz



