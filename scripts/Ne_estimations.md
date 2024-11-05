# Estimate Ne

### currentNe

use only the major chromosomes, so I know the actual distance between snps.

~/spermWhaleRad/analysis/Ne/freebayes_all_chr.vcf.gz
zcat ~/spermWhaleRad/analysis/Ne/freebayes_all_chr.vcf.gz > ~/spermWhaleRad/analysis/Ne/freebayes_all_chr.vcf

subset the vcf to each of the 4 pops:

check this file: `variants_all_chr_NeFormat.ped` bc different results in current Ne than vcf despite same number of snps, etc.

```bash

module load bio/plink/1.90b6.23 

cd ~/spermWhaleRad/analysis/Ne

# pull out only the assembled chromosomes and not scaffolds
# all assembled chr start with NC_
zcat ~/spermWhaleRad/analysis/freebayes/filtered.final.vcf.gz | grep '^#' | grep -v "##contig=<ID=NW_" > ~/spermWhaleRad/analysis/freebayes/vcf_header.txt
zcat ~/spermWhaleRad/analysis/freebayes/filtered.final.vcf.gz | grep '^NC_' | grep -v 'NC_002503.2' > ~/spermWhaleRad/analysis/freebayes/vcf_chrOnly.txt
cat ~/spermWhaleRad/analysis/freebayes/vcf_header.txt ~/spermWhaleRad/analysis/freebayes/vcf_chrOnly.txt > ~/spermWhaleRad/analysis/freebayes/filtered.chrOnly.vcf


~/spermWhaleRad/analysis/Ne/currentNe/currentNe ~/spermWhaleRad/analysis/freebayes/filtered.chrOnly.vcf 21 -o ~/spermWhaleRad/analysis/Ne/chr_only_noK.out

# 4 sibling pairs. 73 total indivs. 
# This means 8 indivs each have value of 1. all other 0. mean(c(rep(1, 8), rep(0, 65))) = 0.11

~/spermWhaleRad/analysis/Ne/currentNe/currentNe ~/spermWhaleRad/analysis/freebayes/filtered.chrOnly.vcf 21 -k -0.11 -o ~/spermWhaleRad/analysis/Ne/chr_only_withK.out


```


# troubleshooting the bug when converting to PED


grab first 5 chrmoosomes:

```bash
zcat ~/spermWhaleRad/analysis/freebayes/filtered.final.vcf.gz | grep -v "^#" |  grep -E "NC_041214|NC_041215|NC_041216|NC_041217|NC_041218"  > ~/spermWhaleRad/analysis/freebayes/subset_test.txt

cat ~/spermWhaleRad/analysis/freebayes/vcf_header.txt ~/spermWhaleRad/analysis/freebayes/subset_test.txt  > ~/spermWhaleRad/analysis/freebayes/subset_test.vcf

sed -i 's/NC_041214.2/1/g' ~/spermWhaleRad/analysis/freebayes/subset_test.vcf
sed -i 's/NC_041215.1/2/g' ~/spermWhaleRad/analysis/freebayes/subset_test.vcf
sed -i 's/NC_041216.1/3/g' ~/spermWhaleRad/analysis/freebayes/subset_test.vcf
sed -i 's/NC_041217.1/4/g' ~/spermWhaleRad/analysis/freebayes/subset_test.vcf
sed -i 's/NC_041218.1/5/g' ~/spermWhaleRad/analysis/freebayes/subset_test.vcf

# convert to plink
module load bio/bcftools/1.11
module load bio/plink
module load lib64/htslib/1.11

cat ~/spermWhaleRad/analysis/freebayes/subset_test.vcf | bgzip  > ~/spermWhaleRad/analysis/freebayes/subset_test.vcf.gz
tabix -f ~/spermWhaleRad/analysis/freebayes/subset_test.vcf.gz
bcftools annotate ~/spermWhaleRad/analysis/freebayes/subset_test.vcf.gz -x ID -I +'%CHROM:%POS' >  ~/spermWhaleRad/analysis/Ne/subset_test_ids.vcf

# conversion with plink
plink --vcf  ~/spermWhaleRad/analysis/Ne/subset_test_ids.vcf --allow-extra-chr --recode --out plink1_conversion

~/plink2 --vcf  ~/spermWhaleRad/analysis/Ne/subset_test_ids.vcf --export ped --out plink2_conversion

# conversion with vcftools
module load bio/vcftools/0.1.16
vcftools --vcf ~/spermWhaleRad/analysis/Ne/subset_test_ids.vcf --out vcftools_conversion --plink


# run currentNe on all
~/spermWhaleRad/analysis/Ne/currentNe/currentNe subset_test_ids.vcf 5 -o  vcf_raw.out
~/spermWhaleRad/analysis/Ne/currentNe/currentNe vcftools_conversion.ped 5 -o  vcftools_conversion.out
~/spermWhaleRad/analysis/Ne/currentNe/currentNe plink2_conversion.ped 5 -o plink2_conversion.out
~/spermWhaleRad/analysis/Ne/currentNe/currentNe plink1_conversion.ped 5 -o plink1_conversion.out


```