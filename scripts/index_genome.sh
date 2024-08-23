#!/bin/bash
#SBATCH -D /home/rbrennan/spermWhaleRad/logout
#SBATCH --mail-type=END
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --partition=standard
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=80G
#SBATCH --time=01:00:00
#SBATCH --job-name=fasta_index
#SBATCH --output=%x.%A.%a.out
#SBATCH --error=%x.%A.%a.err

module load bio/samtools/1.19

cd ~/reference_genomes/sperm_whale/

samtools faidx GCF_002837175.3_ASM283717v5_genomic.fna
