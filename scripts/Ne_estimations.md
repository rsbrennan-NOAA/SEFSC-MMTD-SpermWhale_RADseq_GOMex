# Estimate Ne

## NeEstimator

```r
#Neestimator

# using dartR
library(dartR)

genl <- gl.read.vcf("filtered.final.vcf.gz")
library(dartRverse)

pops <- read.csv("../SW_Metadata.csv")

genl@ind.names <- gsub("b", "",genl@ind.names)

ids <- data.frame(IDs = genl@ind.names)
result <- ids %>%
  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location), by = c("IDs" = "Lab.ID.."))

#genl@pop <- as.factor(result$Pop.Structure.Location)

neest.path.binaries = "C:/Users/Reid.Brennan/Documents/NEestimator2/"

pops <- genl
genl@pop
genl@loc.names <- gsub("\\.", "", genl@loc.names)


nes <- dartR.popgen::gl.LDNe(genl,  outfile = "LD_Ne.txt", neest.path = neest.path.binaries,
               pairing="separate",
               critical = c(0, 0.05), singleton.rm = TRUE, mating = "random")


```



### currentNe

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



# gulf only
grep -v "Atlantic" ~/spermWhaleRad/scripts/populations.txt |  grep -v "Tortuga"   | cut -f 1 > ~/spermWhaleRad/scripts/gulf_only.txt

~/spermWhaleRad/analysis/freebayes/filtered.chrOnly.vcf

vcftools --vcf ~/spermWhaleRad/analysis/freebayes/filtered.chrOnly.vcf \
         --keep ~/spermWhaleRad/scripts/gulf_only.txt \
         --recode \
         --recode-INFO-all \
         --out ~/spermWhaleRad/analysis/freebayes/filtered.chrOnly.gulf_only

~/spermWhaleRad/analysis/Ne/currentNe_new/currentNe ~/spermWhaleRad/analysis/freebayes/filtered.chrOnly.gulf_only.recode.vcf 21 -o ~/spermWhaleRad/analysis/Ne/gulf_only_noK.out

# 41 indivs, mean(c(rep(1, 8), rep(0, 33))) = 0.11
~/spermWhaleRad/analysis/Ne/currentNe_new/currentNe_new ~/spermWhaleRad/analysis/freebayes/filtered.chrOnly.gulf_only.recode.vcf 21 -k -0.195 -o ~/spermWhaleRad/analysis/Ne/gulf_only_withK.out

# atlantic only

grep -v "GOMex" ~/spermWhaleRad/scripts/populations.txt  | cut -f 1 > ~/spermWhaleRad/scripts/atl_only.txt

vcftools --vcf ~/spermWhaleRad/analysis/freebayes/filtered.chrOnly.vcf \
         --keep ~/spermWhaleRad/scripts/atl_only.txt \
         --recode \
         --recode-INFO-all \
         --out ~/spermWhaleRad/analysis/freebayes/filtered.chrOnly.atl_only

~/spermWhaleRad/analysis/Ne/currentNe_new/currentNe ~/spermWhaleRad/analysis/freebayes/filtered.chrOnly.atl_only.recode.vcf 21 -o ~/spermWhaleRad/analysis/Ne/atl_only_noK.out

# 41 indivs, mean(c(rep(1, 8), rep(0, 33))) = 0.11
~/spermWhaleRad/analysis/Ne/currentNe_new/currentNe_new ~/spermWhaleRad/analysis/freebayes/filtered.chrOnly.gulf_only.recode.vcf 21 -k -0.195 -o ~/spermWhaleRad/analysis/Ne/gulf_only_withK.out


```
