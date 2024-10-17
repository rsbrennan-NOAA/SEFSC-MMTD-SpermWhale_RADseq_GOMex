# diversity via pixy
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

setwd("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis/diversity/")

dat <- read.csv("50kb_pi.txt", header=T, sep="\t")
dat.na <- dat[!is.na(dat$avg_pi),]
nrow(dat.na)
# 51748

# use only the assembled chromosomes. 

dat.sub <- dat.na[grep("NC_",dat.na$chromosome),]

nrow(dat.sub)
table(dat.sub$pop)
# 11675 for each pop

hist(dat.sub$no_sites)

ggplot(dat.sub, aes(x=pop, y=(avg_pi))) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 0.002))
sum(dat.na$avg_pi == 0)
# 22275

# find average:
# Note that if the user wishes to combine information across windows 
# (e.g. by averaging) after the fact, they should sum the raw counts 
# and recompute the differences/comparisons ratios, 
# and not take an average of the summary statistics themselves.

sum(dat.sub$count_diffs)/sum(dat.sub$count_comparisons)

mean(dat$avg_pi, na.rm=T)
# actually pretty similar.

pi <- dat.sub %>%
  group_by(pop) %>%
  summarise(overall_pi = sum(count_diffs) / sum(count_comparisons))
pi
#  pop         overall_pi
# Atlantic       0.00137
# Dry Tortuga    0.00133
# NGOMex         0.00140
# WGOMex         0.00132

# two-sided Mann–Whitney–Wilcoxon test with Bonferroni correction

#https://genome.cshlp.org/content/32/10/1952.full.pdf



#####
# bootstrap CI

# Load required libraries
library(dplyr)
library(boot)
library(purrr)

# function to calculate pi
calculate_pi <- function(data, indices) {
  d <- data[indices,] # need to specify indices here bc this is how boot samples. 
  pi <- sum(d$count_diffs) / sum(d$count_comparisons)
  return(pi)
}

# Function to run bootstrap 
bootstrap_population <- function(pop_data) {
  set.seed(123) 
  boot_results <- boot(data = pop_data, statistic = calculate_pi, R = 10000)
  ci <- boot.ci(boot_results, type = "perc")
  
  return(list(
    population = unique(pop_data$pop),
    estimate = boot_results$t0,
    ci_lower = ci$percent[4],
    ci_upper = ci$percent[5]
  ))
}

# run bootstrap
results <- dat.sub %>%
  group_by(pop) %>%
  group_split() %>% # create of df for each pop
  map_dfr(bootstrap_population) # apply function to each df in list


results


#population  estimate ci_lower ci_upper
#<chr>          <dbl>    <dbl>    <dbl>
# Atlantic     0.00137  0.00133  0.00141
# Dry Tortuga  0.00133  0.00130  0.00137
# NGOMex       0.00140  0.00136  0.00143
# WGOMex       0.00132  0.00128  0.00136


# plot

piplot <- ggplot(results, aes(x=population, y=estimate, fill=population)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  geom_point(size=4, pch=21) +
  theme_classic(base_size = 14) +
  ylab("Tajima's pi") +
  scale_fill_manual(values=c("#eac435","#557fc3", "#03cea4", "#fb4d3d"))+
  xlab("")  +
  theme(axis.text.x = element_text(angle = 45, hjust=1))  +
  theme(legend.title=element_blank()) 

ggsave("../figures/pi.png", h=5, w=6)





###### -------------------------------------------------------------------------
# re-run with 10kb windows

dat <- read.csv("10kb_pi.txt", header=T, sep="\t")
dat.na <- dat[!is.na(dat$avg_pi),]
nrow(dat.na)
# 63036

# use only the assembled chromosomes. 

dat.sub <- dat.na[grep("NC_",dat.na$chromosome),]

nrow(dat.sub)
table(dat.sub$pop)
# 14107 for each pop

hist(dat.sub$no_sites)

ggplot(dat.sub, aes(x=pop, y=(avg_pi))) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 0.002))
sum(dat.na$avg_pi == 0)
# 30259

# find average:
# Note that if the user wishes to combine information across windows 
# (e.g. by averaging) after the fact, they should sum the raw counts 
# and recompute the differences/comparisons ratios, 
# and not take an average of the summary statistics themselves.

sum(dat.sub$count_diffs)/sum(dat.sub$count_comparisons)
# 0.00136
mean(dat$avg_pi, na.rm=T)
# 0.00138 

pi <- dat.sub %>%
  group_by(pop) %>%
  summarise(overall_pi = sum(count_diffs) / sum(count_comparisons))
pi
#  pop         overall_pi
# Atlantic       0.00137
# Dry Tortuga    0.00133
# NGOMex         0.00140
# WGOMex         0.00132


#####
# bootstrap CI

# Load required libraries
library(dplyr)
library(boot)
library(purrr)

# function to calculate pi
calculate_pi <- function(data, indices) {
  d <- data[indices,] # need to specify indices here bc this is how boot samples. 
  pi <- sum(d$count_diffs) / sum(d$count_comparisons)
  return(pi)
}

# Function to run bootstrap 
bootstrap_population <- function(pop_data) {
  set.seed(123) 
  boot_results <- boot(data = pop_data, statistic = calculate_pi, R = 10000)
  ci <- boot.ci(boot_results, type = "perc")
  
  return(list(
    population = unique(pop_data$pop),
    estimate = boot_results$t0,
    ci_lower = ci$percent[4],
    ci_upper = ci$percent[5]
  ))
}

# run bootstrap
results <- dat.sub %>%
  group_by(pop) %>%
  group_split() %>% # create of df for each pop
  map_dfr(bootstrap_population) # apply function to each df in list


results


#population  estimate ci_lower ci_upper
#<chr>          <dbl>    <dbl>    <dbl>
# Atlantic     0.00137  0.00134  0.00141
# Dry Tortuga  0.00133  0.00130  0.00137
# NGOMex       0.00140  0.00136  0.00143
# WGOMex       0.00132  0.00128  0.00136



#----------------------------------------------------------------------------------
# obs and het, from stacks
setwd("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis/")

dat <- read.csv("freebayes_filtered_nonrelated.p.sumstats_summary.tsv", header=T, sep="\t")
head(dat)
dat$population <- dat$PopID

obshetplot <- ggplot(dat, aes(x=PopID, y=Obs_Het, fill=population)) +
  geom_errorbar(aes(ymin = Obs_Het-StdErr.2 , ymax = Obs_Het+StdErr.2 ), width = 0.2) +
  geom_point(size=4, pch=21) +
  theme_classic(base_size = 14) +
  ylab("Observed Heterozygosity") +
  xlab("") +
  scale_fill_manual(values=c("#eac435","#557fc3", "#03cea4", "#fb4d3d")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))   +
  theme(legend.title=element_blank())
  

ggsave("../figures/het_obs.png", h=5, w=5)


# expected:

exphetplot <- ggplot(dat, aes(x=PopID, y=Exp_Het, fill=population )) +
  geom_errorbar(aes(ymin = Exp_Het-StdErr.4 , ymax = Exp_Het+StdErr.4 ), width = 0.2) +
  geom_point(size=4, pch=21) +
  theme_classic(base_size = 14) +
  ylab("Expected Heterozygosity") +
  xlab("") +
  scale_fill_manual(values=c("#eac435","#557fc3", "#03cea4", "#fb4d3d")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))  +
  theme(legend.title=element_blank())
ggsave("../figures/het_exp.png", h=5, w=5)


ggarrange(piplot, obshetplot, exphetplot , nrow=1, common.legend =T, legend="top")

ggsave("../figures/diversity.png", h=5, w=10)


# -----------------------------------------------------------------------------------
# check that estimates are consistent with other methods

library(dartR)

genl <- gl.read.vcf("filtered.final.vcf.gz")

pops <- read.csv("../SW_Metadata.csv")

genl@ind.names <- gsub("b", "",genl@ind.names)

ids <- data.frame(IDs = genl@ind.names)
result <- ids %>%
  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location), by = c("IDs" = "Lab.ID.."))

genl@pop <- as.factor(result$Pop.Structure.Location)

gl.report.heterozygosity(genl, method="pop") 
gl.test.heterozygosity(genl, nreps=1000)


# -------------------------------------------------------------------------------
# Fst
# -------------------------------------------------------------------------------

library(snpR)
setwd("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis/")
dat <- snpR::read_vcf("filtered.final.vcf.gz")
pops <- read.csv("../SW_Metadata.csv")

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
my.dat <- calc_pairwise_fst(dat, facets="pop", method = "WC")

# first look at the results, just the head
head(get.snpR.stats(my.dat, facets = "pop", stats = "fst"))


snpout <- get.snpR.stats(my.dat, facets = "pop", stats = "fst")
snpout$fst.matrix$pop
mean(filter(snpout$pairwise, comparison=="Atlantic~Dry_Tortuga")$fst, na.rm=T)
tmp <- filter(snpout$pairwise, comparison=="Atlantic~Dry_Tortuga")$fst
tmp[which(tmp <0)] <- 0
mean(tmp, na.rm=T)

# heatmap of the fst estimates:
png("../figures/fst_heatmap.png", h=4, w=4, units="in", res=300)
plot_pairwise_fst_heatmap(my.dat, facets="pop")
dev.off()


dat <- calc_ho(dat, facets="pop")
dat <- calc_fis(dat, facets="pop")

get.snpR.stats(dat, facets = "pop", stats="ho")
get.snpR.stats(dat, facets = "pop", stats="fis")
#facet    subfacet snp.facet snp.subfacet weighted_mean_fis
#1   pop    Atlantic     .base        .base        0.10168761
#2   pop Dry_Tortuga     .base        .base        0.06291019
#3   pop      NGOMex     .base        .base        0.06533394
#4   pop      WGOMex     .base        .base        0.05719921


#---------------------------------------------
# fst by sexes:meta()
meta(dat)

snp.meta(dat)$CHROM <- gsub("\\.1","",snp.meta(dat)$CHROM)
snp.meta(dat)$CHROM <- gsub("\\.2","",snp.meta(dat)$CHROM)

# calculate fst between the populations
my.dat <- calc_pairwise_fst(dat, facets="sex", method = "WC")

# first look at the results, just the head
head(get.snpR.stats(my.dat, facets = "sex", stats = "fst"))


snpout <- get.snpR.stats(my.dat, facets = "sex.CHROM", stats = "fst")$pairwise
table(snpout$CHROM)
snpout <- snpout[grep("NC_", snpout$CHROM),]
chr1 <- calc_smoothed_averages(x = my.dat, 
                               facets = "sex.CHROM",
                               sigma = 50, # using a window size of 50 kb
                               step = 10) # using a step size of 10kb between windows

# pull out the smoothed values
fst_smooth <- get.snpR.stats(chr1, facets = "sex.CHROM", stats ="fst")$pairwise.window
fst_smooth <- fst_smooth[grep("NC_", fst_smooth$snp.subfacet),]

# make plot
p <- ggplot(snpout, aes(x = position, y = fst, colour = CHROM)) + 
  #geom_point(alpha = 0.5) + 
  theme_bw() +
  theme(legend.position="none") +
  ylim(-0.04,1)+  
  geom_line(data=fst_smooth,
              aes(x=position,y=fst), color="black") +
  facet_wrap(~snp.subfacet, scales = "free_x")

p






chrplot <- snpout$pairwise[grep("NC",snpout$pairwise$CHROM),]
plot_manhattan(chrplot, "fst", chr="CHROM")

calc_smoothed_averages()

sub7 <- calc_smoothed_averages(x = sub7, facets = "pop.chr",
                               sigma = 200, # using a window size of 200
                               step = 50) # using a step size of 50 between windows

snpout$fst.matrix$pop
mean(filter(snpout$pairwise, comparison=="Atlantic~Dry_Tortuga")$fst, na.rm=T)
tmp <- filter(snpout$pairwise, comparison=="Atlantic~Dry_Tortuga")$fst
tmp[which(tmp <0)] <- 0
mean(tmp, na.rm=T)

# heatmap of the fst estimates:
png("../figures/fst_heatmap.png", h=4, w=4, units="in", res=300)
plot_pairwise_fst_heatmap(my.dat, facets="pop")
dev.off()


dat <- calc_ho(dat, facets="pop")
dat <- calc_fis(dat, facets="pop")

get.snpR.stats(dat, facets = "pop", stats="ho")
get.snpR.stats(dat, facets = "pop", stats="fis")
#facet    subfacet snp.facet snp.subfacet weighted_mean_fis
#1   pop    Atlantic     .base        .base        0.10168761
#2   pop Dry_Tortuga     .base        .base        0.06291019
#3   pop      NGOMex     .base        .base        0.06533394
#4   pop      WGOMex     .base        .base        0.05719921



