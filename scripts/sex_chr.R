
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# fst to find sex chromosome.

library(snpR)
library(ggplot2)
library(dplyr)
dat <- snpR::read_vcf("analysis/filtered.final.vcf.gz")
pops <- read.csv("SW_Metadata.csv")

# add meta data information:
## population
colnames(dat) <- gsub("b", "",colnames(dat))
ids <- data.frame(IDs = colnames(dat))
result <- ids %>%
  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location, Sex), by = c("IDs" = "Lab.ID.."))
result$Pop.Structure.Location[result$Pop.Structure.Location =="Dry Tortuga"] <- c("Dry_Tortuga")

sample_meta <- data.frame(pop = result$Pop.Structure.Location, sex=result$Sex)
## order the population
#sample_meta$pop <- factor(sample_meta$pop, levels=c("GA", "HP", "BC", "PC", "TR")) 

# assign meta data to dat
sample.meta(dat) <- sample_meta

# calculate fst between the populations
my.dat <- calc_pairwise_fst(dat, facets="sex", method = "WC")
fst_pvals <- get.snpR.stats(my.dat, facets = "sex", stats = "fst")

fst_pvals$pairwise <- fst_pvals$pairwise %>%
  filter(grepl("^NC", CHROM)) %>%
  mutate(CHROM = as.factor(CHROM)) %>%
  arrange(CHROM, position)

p <- ggplot(fst_pvals$pairwise, aes(x=seq_along(fst), y=fst, color=as.factor(CHROM))) +
  geom_point()+
  scale_color_manual(values = rep(c("dodgerblue4", "black"), length.out = nlevels(fst_pvals$pairwise$CHROM))) +
  theme(legend.position = "none") +
  xlab("SNP position")

ggsave(filename = "figures/male_vs_female_fst.png", p, h=4, w=7)

table(fst_pvals$pairwise$CHROM)


