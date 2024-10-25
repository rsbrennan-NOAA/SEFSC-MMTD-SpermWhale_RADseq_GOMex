#!/bin/bash
#SBATCH --job-name=pixy
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/spermWhaleRad/logout
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 14
#SBATCH --mem=60G
#SBATCH --partition=standard
#SBATCH --time=20:00:00

source ~/.bashrc

mamba activate pixy

cd ~/spermWhaleRad/analysis/diversity/


pixy --stats pi \
--vcf ../freebayes/combined_filtered_invariant.vcf.gz \
--populations ~/spermWhaleRad/scripts/populations.txt \
--window_size 1000 \
--n_cores 14 \
--output_prefix 1kb

# run as one population, bc not really justified to split as there isn't pop structure

awk 'BEGIN{FS=OFS="\t"} {$2="all"; print}' ~/spermWhaleRad/scripts/populations.txt > ~/spermWhaleRad/scripts/populations_all.txt

pixy --stats pi fst dxy \
--vcf ../freebayes/combined_filtered_invariant.vcf.gz \
--populations ~/spermWhaleRad/scripts/populations_all.txt \
--window_size 1000 \
--n_cores 14 \
--output_prefix 1kb_all