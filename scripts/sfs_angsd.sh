#!/bin/bash
#SBATCH --job-name=admixture
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/spermWhaleRad/logout
#SBATCH -c 12
#SBATCH --partition=standard
#SBATCH --time=24:00:00


module load bio/angsd/0.940
#angsd -vcf-gl ~/spermWhaleRad/analysis/freebayes/filtered.final.noMAF.vcf.gz -doSaf 1 -anc /home/rbrennan//reference_genomes/sperm_whale/GCF_002837175.3_ASM283717v5_genomic.fna -out folded

# following https://github.com/z0on/2bRAD_denovo/blob/master/2bRAD_README.sh
cd ~/spermWhaleRad/analysis/demographic_history

FILTERS='-uniqueOnly 1 -skipTriallelic 1 -minMapQ 30 -minQ 30 -maxHetFreq 0.5'
GENOME_REF=/home/rbrennan/reference_genomes/sperm_whale/GCF_002837175.3_ASM283717v5_genomic.fna
TODO="-doHWE 1 -doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 2 -dosnpstat 1 -doGeno 11 -doGlf 2 -anc $GENOME_REF -ref $GENOME_REF"

#angsd -b ~/spermWhaleRad/analysis/bamlist.txt -GL 1 -p 12 -minInd 60 $FILTERS $TODO -underFlowProtect 1 -out folded

realSFS folded.saf.idx  -bootstrap 5 -fold 1 -P 1 >folded.boots

#realSFS folded.saf.idx -ref $GENOME_REF -anc $GENOME_REF -bootstrap 5 -P 1 -fold 1 -resample_chr 1 > folded.boots

# averaging bootstraps, writing sfs file
cat folded.boots | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' > folded.sfs

