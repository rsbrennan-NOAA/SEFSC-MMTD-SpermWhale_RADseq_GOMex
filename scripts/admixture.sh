#!/bin/bash
#SBATCH --job-name=admixture
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/spermWhaleRad/logout
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --time=01:00:00


module load  bio/admixture/1.3.0
module load bio/vcftools/0.1.16
module load bio/plink/1.90b6.23

# need to fix chromosome ids to work with admixture. very annoying. but they need to be numbers
cd ~/spermWhaleRad/analysis/pop_structure

zcat ~/spermWhaleRad/analysis/freebayes/freebayes_ldthin_unrelated.vcf.gz| \
grep -v "^##contig=<ID=N[CW]_"  | awk -F'\t' 'BEGIN {OFS="\t"} 
    /^#/ {print; next} 
    {split($1, a, /\./); sub(/^N[CW]_/, "", a[1]); $1 = a[1]; print}' > LDthin_unrelated_numCorrect.vcf

plink --vcf LDthin_unrelated_numCorrect.vcf \
--make-bed --out LDthin_unrelated_numCorrect \
--allow-extra-chr --double-id


for K in 1 2 3 4 5 6 7 8; \
do admixture --cv LDthin_unrelated_numCorrect.bed $K | tee log${K}.out; done

grep -h CV log*.out | cut -f 3- -d " " > cv.txt



# only females:

vcftools --vcf LDthin_numCorrect_unrelated.vcf --remove ~/spermWhaleRad/analysis/males_only.txt  --recode --recode-INFO-all --stdout  > LDthin_numCorrect_females.vcf

plink --vcf LDthin_numCorrect_females.vcf \
--make-bed --out LDthin_numCorrect_females \
--allow-extra-chr --double-id


for K in 1 2 3 4 5 6 7 8; \
do admixture --cv LDthin_numCorrect_females.bed $K | tee log${K}.out; done

grep -h CV log*.out | cut -f 3- -d " " > cv_females.txt
