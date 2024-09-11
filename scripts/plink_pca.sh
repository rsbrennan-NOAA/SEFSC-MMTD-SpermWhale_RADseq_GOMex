
zcat ~/spermWhaleRad/analysis/freebayes/filtered.final.vcf.gz | grep -v '^#' | cut -f 1 > vcf.1

zcat ~/spermWhaleRad/analysis/freebayes/filtered.final.vcf.gz | grep -v '^#' | cut -f 2 > vcf.2

paste vcf.1 vcf.2 -d ":" > snp_ids.txt

zcat ~/spermWhaleRad/analysis/freebayes/filtered.final.vcf.gz | grep -v '^#' | cut -f 1-2 > vcf.1

zcat ~/spermWhaleRad/analysis/freebayes/filtered.final.vcf.gz | grep -v '^#' | cut -f 4- > vcf.2

paste vcf.1 snp_ids.txt | paste - vcf.2 > vcf.out 

# add header
zcat ~/spermWhaleRad/analysis/freebayes/filtered.final.vcf.gz | grep '^#' | cat - vcf.out | bgzip > ~/spermWhaleRad/analysis/freebayes/filtered.final_ids.vcf.gz

module load bio/plink/1.90b6.23 

plink --vcf ~/spermWhaleRad/analysis/freebayes/filtered.final_ids.vcf.gz \
--indep-pairwise 50 5 0.2 --allow-extra-chr --double-id \
--out variants_pruned

# make thinned vcf

module load bio/vcftools/0.1.16

vcftools --gzvcf ~/spermWhaleRad/analysis/freebayes/filtered.final_ids.vcf.gz --snps variants_pruned.prune.in  --recode --recode-INFO-all --stdout |  bgzip > ~/spermWhaleRad/analysis/freebayes/filtered.final_LDthin.vcf.gz


# make plink files

plink --vcf ~/spermWhaleRad/analysis/freebayes/filtered.final.vcf.gz \
--extract variants_pruned.prune.in \
--make-bed --out variants_NoLD \
--allow-extra-chr --double-id

plink --bfile variants_NoLD \
  --pca --out variants_NoLD_PCA --allow-extra-chr --double-id