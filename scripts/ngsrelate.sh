
module load bio/ngsrelate/20210126

# -c 1 I already called genotypes. 
# -T use genotypes not likelihoods

cd ~/spermWhaleRad/analysis/relatedness

ngsrelate -h ~/spermWhaleRad/analysis/freebayes/filtered.final_LDthin.vcf.gz -T GT -c 1 -O relatedness.res

zcat ~/spermWhaleRad/analysis/freebayes/filtered.final_LDthin.vcf.gz | grep -v ^## | head -n 1 | cut -f 10- | sed 's/\t\t*/\n/g' > indivs.txt
