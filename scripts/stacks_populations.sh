#!/bin/bash
#SBATCH --job-name=populations
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/spermWhaleRad/logout
#SBATCH -c 14
#SBATCH --mem=20G
#SBATCH --partition=standard
#SBATCH --time=2:00:00


module load bio/stacks/2.65

populations --in-vcf ~/spermWhaleRad/analysis/freebayes/filtered.final.vcf.gz --popmap ~/spermWhaleRad/scripts/populations.txt --fstats --threads 14 --out-path ~/spermWhaleRad/analysis/pop_structure/populations_pops
populations --in-vcf ~/spermWhaleRad/analysis/freebayes/filtered.final.vcf.gz --popmap ~/spermWhaleRad/scripts/populations_all.txt --fstats --threads 14 --out-path ~/spermWhaleRad/analysis/pop_structure/populations_all




