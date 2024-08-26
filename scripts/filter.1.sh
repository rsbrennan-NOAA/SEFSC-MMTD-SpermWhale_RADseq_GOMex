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
#SBATCH --time=2:00:00

source ~/.bashrc

mamba activate vcflib-1.0.9

module load bio/vcftools/0.1.16

vcfallelicprimitives ~/spermWhaleRad/analysis/freebayes/freebayes.vcf.gz --keep-info --keep-geno | \
	vcfstreamsort > ~/spermWhaleRad/analysis/freebayes/filtered.1.vcf

# these could be piped, but I want to see the effects of the previous step. not a big deal here bc the vcf is small

vcftools --vcf ~/spermWhaleRad/analysis/freebayes/filtered.1.vcf \
        --max-missing 0.4 --mac 3 --remove-indels --max-alleles 2 --min-alleles 2 --minQ 30  \
        --recode \
        --stdout | bgzip > ~/spermWhaleRad/analysis/freebayes/filtered.2.vcf.gz
