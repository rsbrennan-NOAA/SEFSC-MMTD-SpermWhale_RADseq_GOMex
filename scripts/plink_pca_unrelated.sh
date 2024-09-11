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

# ld thin

# make plink files

plink --vcf ~/spermWhaleRad/analysis/freebayes/freebayes_ldthin_unrelated.vcf.gz \
--make-bed --out variants_unrelated_NoLD \
--allow-extra-chr --double-id

plink --bfile variants_unrelated_NoLD \
  --pca --out variants_NoLD_unrelated_PCA --allow-extra-chr --double-id
