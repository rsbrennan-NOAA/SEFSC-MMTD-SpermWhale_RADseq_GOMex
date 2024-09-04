#!/bin/bash
#SBATCH --job-name=freebayes_allSites
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/spermWhaleRad/logout
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=90G
#SBATCH --partition=standard
#SBATCH --time=72:00:00

##################################
# run freebayes to call variants with ALL SITES!!!!
##################################

# note that this should result in some large files.

#load software---------------------------------------------------------------
source ~/.bashrc

mamba activate freebayes-1.3.6
#module load bio/samtools

#cd ~/spermWhaleRad/analysis/aligned

#samtools index *.bam

#input, output files, directories--------------------------------------------
INDIR=~/spermWhaleRad/analysis/aligned
OUTDIR=~/spermWhaleRad/analysis/freebayes

#reference genome
GEN=~/reference_genomes/sperm_whale/GCF_002837175.3_ASM283717v5_genomic.fna

# make a list of bam files
# POPMAP=../meta/popmap.txt

BAMLIST=~/spermWhaleRad/scripts/bam.list

# these are freebayes scripts found in the same location as the executable
MAKEREGIONS=~/bin/freebayes-1.3.8/scripts/fasta_generate_regions.py

REG=$OUTDIR/regions.txt
# can comment out line below, bc already ran
#python3 $MAKEREGIONS ${GEN}.fai 5000000 >$REG

# run freebayes-parallel--------------------------------------------------------

fb_parallel=~/bin/freebayes-1.3.8/scripts/freebayes-parallel

bash $fb_parallel \
        $REG 20 \
        -f ${GEN} \
        --bam-list $BAMLIST \
        -m 30 \
        -q 20 \
	--min-coverage 10 \
	--report-monomorphic \
        --skip-coverage 5000 | \
bgzip -c >$OUTDIR/freebayes_allsites.vcf.gz

# index the vcf file
tabix -p vcf $OUTDIR/freebayes_allsites.vcf.gz
