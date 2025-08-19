
# make nexus file for popart:
library(ape)

myseqs <- read.dna("analysis/mitotyping/Pmac_All_Align_complete_CR_truncated_RADSamples_only-DropAmbig.fasta",format="fasta",as.matrix=FALSE)
myseqs.names <- names(myseqs)

myseqs.nex.fn <- "analysis/mitotyping/Pmac_All.nex" # output in Nexus format
write.nexus.data(myseqs,myseqs.nex.fn,interleaved=FALSE,gap="-",
                 missing="?",datablock = FALSE)


#### ----------------------------------
library(pegas)
library(ape)
library(seqinr)
library(tidyverse)
library(dplyr)
library(ggpubr)

data <- read.FASTA("analysis/mitotyping/Pmac_All_Align_complete_CR_truncated_RADSamples_only-DropAmbig.fasta")

length(data)  
hap <- haplotype(data, strict = FALSE, trailingGapsAsN = FALSE)
sapply(hap, function(x) any(!x %in% c("a", "t", "g", "c", "A", "T", "G", "C")))

names <- data.frame(IDs = labels(data))

# add region:

pops <- read.csv("meta_data_submission.csv")

pops[,c(1,2)]

hapInfo <- stack(setNames(attr(hap,"index"),rownames(hap)))

names(hapInfo) <- c("index","haplotype")
head(hapInfo)
hapInfo$names <- names$IDs[hapInfo$index]
head(hapInfo)
names(hapInfo) <- c("index","haplotype", "IDs")

str(hapInfo)


names_with_pops <- merge(hapInfo, 
                         pops, 
                         by.x = "IDs", 
                         by.y = "Lab.ID..")


# Create the frequency table 
table(names_with_pops$Pop.Structure.Location)
#   Atlantic Dry Tortuga      NGOMex      WGOMex 
#      14          17          20          18 
pop_freq <- table(names_with_pops$hap, names_with_pops$Pop.Structure.Location)


net <- haploNet(hap)


pdf("../figures/mito_network.pdf", h=4, w=4)

plot(net, 
     size = attr(net, "freq"),
     pie = pop_freq,
     scale.ratio=20, cex=0.7,
     show.mutation=1, 
     bg=c("#eac435","#557fc3", "#03cea4", "#fb4d3d"))

legend("bottomleft", colnames(pop_freq), 
       col=c("#eac435","#557fc3", "#03cea4", "#fb4d3d"), pch=19, ncol=1, cex=0.7)

dev.off()

# this plot is very ugly. run it in popart instead and use that for ms. Results are very similar. 



#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# calc pop gen 
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

library(adegenet)
library(hierfstat)
dat.gid <- DNAbin2genind(data)

strata(dat.gid) <- data.frame(Pop.Structure.Location = names_with_pops[,c(4)])
setPop(dat.gid) <- ~Pop.Structure.Location
dat.gid

  
dat.gid@tab

out <- genind2hierfstat(dat.gid)
a <- genet.dist(out, diploid=FALSE, method = "Nei87")
a
#            Atlantic NGOMex WGOMex
#NGOMex        0.3068              
#WGOMex        0.4962 0.0165       
#Dry Tortuga   0.0011 0.4520 0.6415

pairwise.neifst(out,diploid=FALSE) # gives same as above
a <- genet.dist(out, diploid=FALSE, method = "WC84")
a

# Calculate Nei's GST (FST)
nei.fst <- mmod::pairwise_Gst_Nei(dat.gid)

# command below is wc87 fst
wc.fst <- pairwise.WCfst(out, diploid=F)
boot.ppfst(dat=out,nboot=1000,quant=c(0.025,0.975),diploid=FALSE)

# boot.ppfst runs bootstrapping. but I actually want to do permutation, mix pop assignments
# this is what snpR definitely does and, I think arlequin both do.
# makes sense to test null of no structure.

# Calculate observed FST
observed_fst <- pairwise.WCfst(out, diploid=F)
#Atlantic     NGOMex     WGOMex Dry Tortuga
#Atlantic             NA 0.31288267 0.51388926 0.002259552
#NGOMex      0.312882669         NA 0.01510005 0.451089959
#WGOMex      0.513889262 0.01510005         NA 0.643979010
#Dry Tortuga 0.002259552 0.45108996 0.64397901          NA

npop <- nrow(observed_fst)
p_values <- matrix(NA, npop, npop)

nperm <- 1000
perm_dist <- array(NA, dim=c(npop, npop, nperm))

# Run permutations
for(i in 1:nperm) {
  out_perm <- out
  out_perm$pop <- sample(out$pop)
  perm_fst <- pairwise.WCfst(out_perm, diploid=F)
  perm_fst[perm_fst < 0] <- 0
  perm_dist[,,i] <- perm_fst
}

observed_fst[observed_fst < 0] <- 0

# calc pvals
for(i in 1:npop) {
  for(j in 1:npop) {
    p_values[i,j] <- sum(perm_dist[i,j,] >= observed_fst[i,j])/nperm
  }
}

p_values


###--------------------------
# plot
###--------------------------
# convert fst to long

melt_fst <- as_tibble(observed_fst, rownames = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")%>%
  mutate(across(c(Var1, Var2), ~str_replace(., "Dry_Tortuga", "Dry Tortuga")))

melt_fst <- subset(melt_fst, !(Var1 != Var2 & is.na(value)))

melt_fst <- melt_fst %>%
  group_by(value) %>%
  filter(!(duplicated(paste0(pmin(Var1, Var2), pmax(Var1, Var2))) & !is.na(value)))

# convert p-values to long:
rownames(p_values) <- row.names(observed_fst)
colnames(p_values) <- colnames(observed_fst)

melt_fst_pval<- as_tibble(p_values, rownames = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")%>%
  mutate(across(c(Var1, Var2), ~str_replace(., "Dry_Tortuga", "Dry Tortuga")))

melt_fst_pval <- melt_fst_pval%>%
  group_by(value) %>%
  filter(!(duplicated(paste0(pmin(Var1, Var2), pmax(Var1, Var2))) & !is.na(value)))


# flip order for pvals, to put on other triangle:
melt_fst_pval2 <- melt_fst_pval[, c("Var2", "Var1", "value")]
colnames(melt_fst_pval2) <- c("Var1", "Var2", "value")

melt_fst$group <- c("fst")
melt_fst_pval2$group <- c("pval")
fst_plot_dat <- rbind(melt_fst, melt_fst_pval2)

fst_plot_dat$group[melt_fst$Var1 == melt_fst$Var2] <- "diagonal"

fst_plot_dat <- fst_plot_dat %>% 
                    distinct()


fst_plot_dat$Var1 <- factor(fst_plot_dat$Var1, c("Atlantic", "Dry Tortuga", "NGOMex", "WGOMex"))
fst_plot_dat$Var2 <- factor(fst_plot_dat$Var2, c("Atlantic", "Dry Tortuga", "NGOMex", "WGOMex"))


fst_plot_dat2 <- fst_plot_dat

# drop p-vals
fst_plot_dat3 <- subset(fst_plot_dat2, !group %in% c( "pval"))
fst_plot_dat4 <- subset(fst_plot_dat2, group %in% c( "pval"))
fst_plot_diags <- subset(fst_plot_dat2, group %in% c( "diagonal"))

# need to flip the order of a couple for plots:
fst_plot_dat3[c(7,9), c("Var1","Var2")] <- fst_plot_dat3[c(7,9), c("Var2","Var1")]

fst_plot_dat4[c(5,6), c("Var1","Var2")] <- fst_plot_dat4[c(5,6), c("Var2","Var1")]

p_fst <- ggplot(fst_plot_dat3, aes(y = Var2, x = Var1, fill = value)) +
  geom_tile(size=0.2, color="black") +
  geom_text(aes(label = round(value, 3)), na.rm = TRUE, size=3) +
  scale_fill_gradient(limits = c(0, 0.76),
                      low = "#FFCCCC", 
                      high = "firebrick3", 
                      na.value = "white") +  
  theme_bw(base_size = 13) +
  labs(x = NULL, y = NULL, fill =  bquote(F[ST])) +
  geom_tile(data = fst_plot_diags, 
            aes(x = Var1, y = Var2),
            fill="grey",
            size=0.2, color="black") +
  geom_tile(data = fst_plot_dat4, 
            aes(x = Var1, y = Var2),
            fill="white",
            size=0.2, color="black") +
  geom_text(data = fst_plot_dat4, 
            aes(label = round(value, 2)),
            size=3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.position = c(0.95, 0.05),
    #legend.justification = c(1, 0),
    legend.background = element_blank(),  # Remove legend background
    legend.margin = margin(0, 0, 0, 0),  # Remove legend margin
    legend.box.margin = margin(0, 0, 0, -10),  # negative value moves legend left
    legend.key.size = unit(0.5, "lines"),  #  reduced legend key size
    legend.title = element_text(size = 9, margin = margin(b = 0)),  # Reduce bottom margin of title
    legend.text = element_blank(),
    #legend.text = element_text(size = 7, margin = margin(l = 1, r = 0)),  # Reduced legend text size
    legend.spacing.y = unit(0.05, "cm"),  # Reduced spacing between legend elements
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_fst

### ----------------------
# amova:
library(poppr)

dat.gc <- as.genclone(dat.gid)
amova.result <- poppr.amova(dat.gc, ~Pop.Structure.Location,method = "pegas", 
                            within = FALSE)
amova.result <- poppr.amova(dat.gc, ~Pop.Structure.Location, 
                            within = FALSE)
amova.result
amova.test <- randtest(amova.result, nrepet = 999) # Test for significance
plot(amova.test)

amova.test


#### -----------------------------------------------------------
#### -----------------------------------------------------------
#### -----------------------------------------------------------
# nucleotide diversity
#### -----------------------------------------------------------
#### -----------------------------------------------------------
#### -----------------------------------------------------------

mydna <- read.dna("analysis/mitotyping/Pmac_All_Align_complete_CR_truncated_RADSamples_only-DropAmbig.fasta",format="fasta")

ss <- sapply(2:69, function(z){
  length(seg.sites(mydna[1:z, ]))
})
plot(2:69, ss, col = "red", xlab = "No of sequences", ylab = "Segregating sites", las = 1)

seg.sites(mydna)
seg.sites(hap)


nuc.div(mydna)
# [1] 0.001800547

hap.div(mydna)
#0.7374254

# get CI from bootstrapping
n_reps <- 1000
n_seqs <- nrow(mydna)
results <- numeric(n_reps)

for(i in 1:n_reps) {
  # Sample individuals with replacement
  ind_sample <- sample(1:n_seqs, n_seqs, replace=TRUE)
  boot_seqs <- mydna[ind_sample,]
  results[i] <- nuc.div((boot_seqs))
}

ci <- quantile(results, c(0.025, 0.975))
print(paste("Mean pi, boots:", round(mean(results), 5)))
print(paste("95% CI:", round(ci[1], 5), "-", round(ci[2], 5)))
#[1] "95% CI: 0.0014 - 0.0019"

mt_all <- data.frame(
  Population = "All Populations",
  n_samples = 73,
  nuc_div = nuc.div(mydna),
  mean_boot = mean(results), # I know this is wrong
  ci_lower = ci[1],
  ci_upper = ci[2],
  row.names=NULL)


##### --------------------------
# also separate by population:

metadata <- names_with_pops
n_reps=1000

# Split data by population
populations <- unique(metadata$Pop.Structure.Location)

# make empty df to hold results
results_df <- data.frame(
  Population = character(),
  n_samples = numeric(),
  nuc_div = numeric(),
  mean_boot = numeric(),
  ci_lower = numeric(),
  ci_upper = numeric(),
  stringsAsFactors = FALSE
)

for(pop in populations) {
  pop_ids <- metadata$IDs[metadata$Pop.Structure.Location == pop]
  pop_dna <- mydna[pop_ids,]
  
  pop_nuc_div <- nuc.div(pop_dna)
  n_seqs <- length(pop_ids)
  boot_results <- numeric(n_reps)
  
  for(i in 1:n_reps) {
    ind_sample <- sample(1:n_seqs, n_seqs, replace=TRUE)
    boot_seqs <- pop_dna[ind_sample,]
    boot_results[i] <- nuc.div(boot_seqs)
  }
  
  ci <- quantile(boot_results, c(0.025, 0.975))
  
  results_df <- rbind(results_df, data.frame(
    Population = pop,
    n_samples = n_seqs,
    nuc_div = pop_nuc_div,
    mean_boot = mean(boot_results),
    ci_lower = ci[1],
    ci_upper = ci[2],
    row.names=NULL
  ))
}


results_combined <- rbind(results_df, mt_all)
population_order <- c("Atlantic", "Dry Tortuga","All Populations", "NGOMex", "WGOMex")
results_combined$Population <- factor(results_combined$Population, levels = population_order)
legend_order <- c("All Populations", "Atlantic", "Dry Tortuga", "NGOMex", "WGOMex")
results_combined$population_legend <- factor(results_combined$Population, levels = legend_order)

divplot <- ggplot(results_combined, 
                 aes(x = Population, 
                     y = nuc_div, 
                     fill = population_legend, shape=population_legend)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  geom_point(aes(size = ifelse(Population == "All Populations", "All", "Population"))) +
  scale_size_manual(values = c("All" = 8, "Population" = 4), guide = "none") +
  theme_classic(base_size = 12) +
  ylab("Genetic\ndiversity") +
  scale_color_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_fill_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,21,22,23,24))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(
    breaks = c(0, 0.001, 0.002),
    labels = function(x) format(x, scientific = FALSE, big.mark = "")
  )


divplot

ggsave("figures/mt_diversity.png", h=5, w=6)
ggsave("figures/mt_diversity.pdf", h=5, w=6)

#drop legend:
divplot_noleg <-divplot + 
  theme(legend.position = "none")


# combine plots:

p_fst2 <- p_fst + theme(plot.margin = margin(t = 0.1,  # Top margin
                                             r = 0.5,  # Right margin
                                             b = 0.2,  # Bottom margin
                                             l = 0.5,
                                             unit="cm")) # Left margin
divplot2 <- divplot_noleg + theme(plot.margin = margin(t = 0.05,  # Top margin
                                               r = 0.5,  # Right margin
                                               b = 0.1,  # Bottom margin
                                               l = 0.5,
                                               unit="cm")) # Left margin
p_out1 <- ggarrange(p_fst2, divplot2, ncol = 1, nrow = 2, labels=c("B", "C"))

ggsave(file="figures/mito_diversity.pdf", w=3, h=4)

p_fst2 <- p_fst + theme(plot.margin = margin(t = 0.15,  # Top margin
                                             r = 0.1,  # Right margin
                                             b = 0.1,  # Bottom margin
                                             l = 0.1,
                                             unit="cm")) # Left margin
divplot2 <- divplot_noleg + theme(plot.margin = margin(t = 0.15,  # Top margin
                                                       r = 0.1,  # Right margin
                                                       b = 0.0,  # Bottom margin
                                                       l = 0.3,
                                                       unit="cm")) # Left margin
p_out1 <- ggarrange(p_fst2, divplot2, ncol = 2, nrow = 1, labels=c("B", "C"))
ggsave(file="figures/mito_diversity-2.pdf", w=5.5, h=2.5)



########################
# now atlantic vs gulf:

metadata$Region <- ifelse(metadata$Pop.Structure.Location %in% c("Atlantic", "Dry Tortuga"), 
                          "Atlantic_DryTortuga", "Gulf_of_Mexico")

# Run analysis with new grouping
populations <- unique(metadata$Region)

for(pop in populations) {
  pop_ids <- metadata$IDs[metadata$Region == pop]
  pop_dna <- mydna[pop_ids,]
  
  pop_nuc_div <- nuc.div(pop_dna)
  n_seqs <- length(pop_ids)
  boot_results <- numeric(n_reps)
  
  for(i in 1:n_reps) {
    ind_sample <- sample(1:n_seqs, n_seqs, replace=TRUE)
    boot_seqs <- pop_dna[ind_sample,]
    boot_results[i] <- nuc.div(boot_seqs)
  }
  
  ci <- quantile(boot_results, c(0.025, 0.975))
  
  results_df <- rbind(results_df, data.frame(
    Population = pop,
    n_samples = n_seqs,
    nuc_div = pop_nuc_div,
    mean_boot = mean(boot_results),
    ci_lower = ci[1],
    ci_upper = ci[2],
    row.names=NULL
  ))
}








#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# separate males and females


# make nexus file for popart:
library(ape)

#myseqs <- read.dna("analysis/mitotyping/Pmac_All_Align_complete_CR_truncated_RADSamples_only_females.fasta",format="fasta",as.matrix=FALSE)
myseqs <- read.dna("analysis/mitotyping/Pmac_All_Align_complete_CR_truncated_RADSamples_only-DropAmbig.fasta",format="fasta",as.matrix=FALSE)
myseqs.names <- names(myseqs)

female_ids <- pops$Lab.ID..[pops$Sex == "F"]
length(female_ids)  # How many females?

# Subset sequences to only females
female_seqs <- myseqs[names(myseqs) %in% female_ids]
length(female_seqs)  # Should match number of female_ids
# 3 short bc removed the heterplasmic

myseqs.nex.fn <- "analysis/mitotyping/Pmac_All_females.nex" # output in Nexus format
write.nexus.data(as.character(female_seqs),myseqs.nex.fn,interleaved=FALSE,gap="-",
                 missing="?",datablock = FALSE)

hap_original <- haplotype(myseqs, strict = FALSE, trailingGapsAsN = FALSE)

# Female subset  
hap_females <- haplotype(female_seqs, strict = FALSE, trailingGapsAsN = FALSE)



#### ----------------------------------
library(pegas)
library(ape)
library(seqinr)
library(tidyverse)
library(dplyr)
library(ggpubr)


hap <- hap_females

names <- data.frame(IDs = labels(female_seqs))

# add region:

pops <- read.csv("meta_data_submission.csv")

pops[,c(1,2)]

hapInfo <- stack(setNames(attr(hap,"index"),rownames(hap)))

names(hapInfo) <- c("index","haplotype")
head(hapInfo)
hapInfo$names <- names$IDs[hapInfo$index]
head(hapInfo)
names(hapInfo) <- c("index","haplotype", "IDs")

str(hapInfo)

names_with_pops <- merge(hapInfo, 
                         pops, 
                         by.x = "IDs", 
                         by.y = "Lab.ID..")


# Create the frequency table 
table(names_with_pops$Pop.Structure.Location)
#Atlantic    Dry Tortuga      NGOMex      WGOMex 
#     8            9          12          16 
pop_freq <- table(names_with_pops$hap, names_with_pops$Pop.Structure.Location)


net <- haploNet(hap)


pdf("../figures/mito_network_females.pdf", h=4, w=4)

plot(net, 
     size = attr(net, "freq"),
     pie = pop_freq,
     scale.ratio=20, cex=0.7,
     show.mutation=1, 
     bg=c("#eac435","#557fc3", "#03cea4", "#fb4d3d"))

legend("bottomleft", colnames(pop_freq), 
       col=c("#eac435","#557fc3", "#03cea4", "#fb4d3d"), pch=19, ncol=1, cex=0.7)

dev.off()

# this plot is very ugly. run it in popart instead and use that for ms. Results are very similar. 



# pop gen ---------------------------------------------------------------

library(adegenet)
library(hierfstat)
dat.gid <- DNAbin2genind(female_seqs)

strata(dat.gid) <- data.frame(Pop.Structure.Location = names_with_pops[,c(4)])
setPop(dat.gid) <- ~Pop.Structure.Location
dat.gid


dat.gid@tab

out <- genind2hierfstat(dat.gid)
a <- genet.dist(out, diploid=FALSE, method = "Nei87")
a
#            Atlantic NGOMex WGOMex
#NGOMex        0.6575                
#WGOMex        0.6922  0.0562        
#Dry Tortuga  -0.0385  0.6337  0.6685

pairwise.neifst(out,diploid=FALSE) # gives same as above
a <- genet.dist(out, diploid=FALSE, method = "WC84")
a

# Calculate Nei's GST (FST)
nei.fst <- mmod::pairwise_Gst_Nei(dat.gid)

# command below is wc87 fst
wc.fst <- pairwise.WCfst(out, diploid=F)
boot.ppfst(dat=out,nboot=1000,quant=c(0.025,0.975),diploid=FALSE)

# boot.ppfst runs bootstrapping. but I actually want to do permutation, mix pop assignments
# this is what snpR definitely does and, I think arlequin both do.
# makes sense to test null of no structure.

# Calculate observed FST
observed_fst <- pairwise.WCfst(out, diploid=F)
#               Atlantic        NGOMex        WGOMex Dry Tortuga
#Atlantic             NA 0.68586387 0.76139188 -0.03846154
#NGOMex       0.68586387         NA 0.06582713  0.65445026
#WGOMex       0.76139188 0.06582713         NA  0.72844037
#Dry Tortuga -0.03846154 0.65445026 0.72844037          NA

npop <- nrow(observed_fst)
p_values <- matrix(NA, npop, npop)

nperm <- 1000
perm_dist <- array(NA, dim=c(npop, npop, nperm))

# Run permutations
for(i in 1:nperm) {
  out_perm <- out
  out_perm$pop <- sample(out$pop)
  perm_fst <- pairwise.WCfst(out_perm, diploid=F)
  perm_fst[perm_fst < 0] <- 0
  perm_dist[,,i] <- perm_fst
}

observed_fst[observed_fst < 0] <- 0

# calc pvals
for(i in 1:npop) {
  for(j in 1:npop) {
    p_values[i,j] <- (sum(perm_dist[i,j,] >= observed_fst[i,j])+1)/(nperm+1)
  }
}

p_values


###--------------------------
# plot
###--------------------------
# convert fst to long

melt_fst <- as_tibble(observed_fst, rownames = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")%>%
  mutate(across(c(Var1, Var2), ~str_replace(., "Dry_Tortuga", "Dry Tortuga")))

melt_fst <- subset(melt_fst, !(Var1 != Var2 & is.na(value)))

melt_fst <- melt_fst %>%
  group_by(value) %>%
  filter(!(duplicated(paste0(pmin(Var1, Var2), pmax(Var1, Var2))) & !is.na(value)))

# convert p-values to long:
rownames(p_values) <- row.names(observed_fst)
colnames(p_values) <- colnames(observed_fst)

melt_fst_pval<- as_tibble(p_values, rownames = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")%>%
  mutate(across(c(Var1, Var2), ~str_replace(., "Dry_Tortuga", "Dry Tortuga")))

melt_fst_pval <- melt_fst_pval%>%
  group_by(value) %>%
  filter(!(duplicated(paste0(pmin(Var1, Var2), pmax(Var1, Var2))) & !is.na(value)))


# flip order for pvals, to put on other triangle:
melt_fst_pval2 <- melt_fst_pval[, c("Var2", "Var1", "value")]
colnames(melt_fst_pval2) <- c("Var1", "Var2", "value")

melt_fst$group <- c("fst")
melt_fst_pval2$group <- c("pval")
fst_plot_dat <- rbind(melt_fst, melt_fst_pval2)

fst_plot_dat$group[melt_fst$Var1 == melt_fst$Var2] <- "diagonal"

fst_plot_dat <- fst_plot_dat %>% 
  distinct()


fst_plot_dat$Var1 <- factor(fst_plot_dat$Var1, c("Atlantic", "Dry Tortuga", "NGOMex", "WGOMex"))
fst_plot_dat$Var2 <- factor(fst_plot_dat$Var2, c("Atlantic", "Dry Tortuga", "NGOMex", "WGOMex"))


fst_plot_dat2 <- fst_plot_dat

# drop p-vals
fst_plot_dat3 <- subset(fst_plot_dat2, !group %in% c( "pval"))
fst_plot_dat4 <- subset(fst_plot_dat2, group %in% c( "pval"))
fst_plot_diags <- subset(fst_plot_dat2, group %in% c( "diagonal"))

# need to flip the order of a couple for plots:
fst_plot_dat3[c(7,9), c("Var1","Var2")] <- fst_plot_dat3[c(7,9), c("Var2","Var1")]

fst_plot_dat4[c(5,6), c("Var1","Var2")] <- fst_plot_dat4[c(5,6), c("Var2","Var1")]

p_fst <- ggplot(fst_plot_dat3, aes(y = Var2, x = Var1, fill = value)) +
  geom_tile(size=0.2, color="black") +
  geom_text(aes(label = round(value, 2)), na.rm = TRUE, size=3) +
  scale_fill_gradient(limits = c(0, 0.77),
                      low = "#FFCCCC", 
                      high = "firebrick3", 
                      na.value = "white") +  
  theme_bw(base_size = 13) +
  labs(x = NULL, y = NULL, fill =  bquote(F[ST])) +
  geom_tile(data = fst_plot_diags, 
            aes(x = Var1, y = Var2),
            fill="grey",
            size=0.2, color="black") +
  geom_tile(data = fst_plot_dat4, 
            aes(x = Var1, y = Var2),
            fill="white",
            size=0.2, color="black") +
  geom_text(data = fst_plot_dat4, 
            aes(label = round(value, 2)),
            size=3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.position = c(0.95, 0.05),
    #legend.justification = c(1, 0),
    legend.background = element_blank(),  # Remove legend background
    legend.margin = margin(0, 0, 0, 0),  # Remove legend margin
    legend.box.margin = margin(0, 0, 0, -10),  # negative value moves legend left
    legend.key.size = unit(0.5, "lines"),  #  reduced legend key size
    legend.title = element_text(size = 9, margin = margin(b = 0)),  # Reduce bottom margin of title
    legend.text = element_blank(),
    #legend.text = element_text(size = 7, margin = margin(l = 1, r = 0)),  # Reduced legend text size
    legend.spacing.y = unit(0.05, "cm"),  # Reduced spacing between legend elements
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p_fst




#### -----------------------------------------------------------
#### -----------------------------------------------------------
#### -----------------------------------------------------------
# nucleotide diversity

mydna <- read.dna("analysis/mitotyping/Pmac_All_Align_complete_CR_truncated_RADSamples_only-DropAmbig.fasta",
                  format="fasta", as.matrix=TRUE)
female_ids <- pops$Lab.ID..[pops$Sex == "F"]
length(female_ids)

# Subset to females only
female_seqs <- mydna[rownames(mydna) %in% female_ids, ]
dim(female_seqs) 
mydna <- female_seqs


seg.sites(mydna)
seg.sites(hap)

nuc.div(mydna)
# [1] 0.001599238

hap.div(mydna)
#0.6747475

# get CI from bootstrapping
n_reps <- 1000
n_seqs <- nrow(mydna)
results <- numeric(n_reps)

for(i in 1:n_reps) {
  # Sample individuals with replacement
  ind_sample <- sample(1:n_seqs, n_seqs, replace=TRUE)
  boot_seqs <- mydna[ind_sample,]
  results[i] <- nuc.div((boot_seqs))
}

ci <- quantile(results, c(0.025, 0.975))
print(paste("Mean pi, boots:", round(mean(results), 5)))
print(paste("95% CI:", round(ci[1], 5), "-", round(ci[2], 5)))
#[1] "95% CI: 0.00114 - 0.00192"

mt_all <- data.frame(
  Population = "All Populations",
  n_samples = 45,
  nuc_div = nuc.div(mydna),
  mean_boot = mean(results), # I know this is wrong
  ci_lower = ci[1],
  ci_upper = ci[2],
  row.names=NULL)


##### --------------------------
# also separate by population:


metadata <- names_with_pops
n_reps=1000

# Split data by population
populations <- unique(metadata$Pop.Structure.Location)

# make empty df to hold results
results_df <- data.frame(
  Population = character(),
  n_samples = numeric(),
  nuc_div = numeric(),
  mean_boot = numeric(),
  ci_lower = numeric(),
  ci_upper = numeric(),
  stringsAsFactors = FALSE
)

for(pop in populations) {
  pop_ids <- metadata$IDs[metadata$Pop.Structure.Location == pop]
  pop_dna <- mydna[pop_ids,]
  
  pop_nuc_div <- nuc.div(pop_dna)
  n_seqs <- length(pop_ids)
  boot_results <- numeric(n_reps)
  
  for(i in 1:n_reps) {
    ind_sample <- sample(1:n_seqs, n_seqs, replace=TRUE)
    boot_seqs <- pop_dna[ind_sample,]
    boot_results[i] <- nuc.div(boot_seqs)
  }
  
  ci <- quantile(boot_results, c(0.025, 0.975))
  
  results_df <- rbind(results_df, data.frame(
    Population = pop,
    n_samples = n_seqs,
    nuc_div = pop_nuc_div,
    mean_boot = mean(boot_results),
    ci_lower = ci[1],
    ci_upper = ci[2],
    row.names=NULL
  ))
}


results_combined <- rbind(results_df, mt_all)
population_order <- c("Atlantic", "Dry Tortuga","All Populations", "NGOMex", "WGOMex")
results_combined$Population <- factor(results_combined$Population, levels = population_order)
legend_order <- c("All Populations", "Atlantic", "Dry Tortuga", "NGOMex", "WGOMex")
results_combined$population_legend <- factor(results_combined$Population, levels = legend_order)

#Population n_samples           nuc_div      mean_boot     ci_lower     ci_upper population_legend
#1        Atlantic         8 0.0015706806 0.0013626776 0.0005974196 0.0017576664          Atlantic
#2          NGOMex        12 0.0004283674 0.0003936855 0.0000000000 0.0005711566            NGOMex
#3          WGOMex        16 0.0001308901 0.0001178534 0.0000000000 0.0003403141            WGOMex
#4     Dry Tortuga         9 0.0015706806 0.0013872019 0.0006980803 0.0017452007       Dry Tortuga
#5 All Populations        45 0.0015992385 0.0015544725 0.0010765773 0.0019250093   All Populations


divplot <- ggplot(results_combined, 
                  aes(x = Population, 
                      y = nuc_div, 
                      fill = population_legend, shape=population_legend)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  geom_point(aes(size = ifelse(Population == "All Populations", "All", "Population"))) +
  scale_size_manual(values = c("All" = 8, "Population" = 4), guide = "none") +
  theme_classic(base_size = 12) +
  ylab("Genetic\ndiversity") +
  scale_color_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_fill_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,21,22,23,24))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(
    breaks = c(0, 0.001, 0.002),
    labels = function(x) format(x, scientific = FALSE, big.mark = "")
  )


divplot

ggsave("figures/mt_diversity_females.png", h=5, w=6)
ggsave("figures/mt_diversity_females.pdf", h=5, w=6)

#drop legend:
divplot_noleg <-divplot + 
  theme(legend.position = "none")


# combine plots:

p_fst2 <- p_fst + theme(plot.margin = margin(t = 0.1,  # Top margin
                                             r = 0.5,  # Right margin
                                             b = 0.2,  # Bottom margin
                                             l = 0.5,
                                             unit="cm")) # Left margin
divplot2 <- divplot_noleg + theme(plot.margin = margin(t = 0.05,  # Top margin
                                                       r = 0.5,  # Right margin
                                                       b = 0.1,  # Bottom margin
                                                       l = 0.5,
                                                       unit="cm")) # Left margin
p_out1 <- ggarrange(p_fst2, divplot2, ncol = 1, nrow = 2, labels=c("B", "C"))

ggsave(file="figures/mito_diversity_females.pdf", w=3, h=4)

p_fst2 <- p_fst + theme(plot.margin = margin(t = 0.15,  # Top margin
                                             r = 0.1,  # Right margin
                                             b = 0.1,  # Bottom margin
                                             l = 0.1,
                                             unit="cm")) # Left margin
divplot2 <- divplot_noleg + theme(plot.margin = margin(t = 0.15,  # Top margin
                                                       r = 0.1,  # Right margin
                                                       b = 0.0,  # Bottom margin
                                                       l = 0.3,
                                                       unit="cm")) # Left margin
p_out1 <- ggarrange(p_fst2, divplot2, ncol = 2, nrow = 1, labels=c("B", "C"))
ggsave(file="figures/mito_diversity-2_females.pdf", w=5.5, h=2.5)












#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Males
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# make nexus file for popart:
library(ape)

myseqs <- read.dna("analysis/mitotyping/Pmac_All_Align_complete_CR_truncated_RADSamples_only-DropAmbig.fasta",format="fasta",as.matrix=FALSE)
myseqs.names <- names(myseqs)

female_ids <- pops$Lab.ID..[pops$Sex == "M"]
length(female_ids)  # How many females?
# 25

# Subset sequences to only females
female_seqs <- myseqs[names(myseqs) %in% female_ids]

female_ids[(!female_ids %in% names(female_seqs))]
length(female_seqs)  # Should match number of female_ids
# but 1 shorter, bc pmac 117 removed bc heteroplasmic. 

myseqs.nex.fn <- "analysis/mitotyping/Pmac_All_males.nex" # output in Nexus format
write.nexus.data(as.character(female_seqs),myseqs.nex.fn,interleaved=FALSE,gap="-",
                 missing="?",datablock = FALSE)

hap_original <- haplotype(myseqs, strict = FALSE, trailingGapsAsN = FALSE)

# Female subset  
hap_females <- haplotype(female_seqs, strict = FALSE, trailingGapsAsN = FALSE)



#### ----------------------------------
library(pegas)
library(ape)
library(seqinr)
library(tidyverse)
library(dplyr)
library(ggpubr)


hap <- hap_females

names <- data.frame(IDs = labels(female_seqs))

# add region:

pops <- read.csv("meta_data_submission.csv")

pops[,c(1,2)]

hapInfo <- stack(setNames(attr(hap,"index"),rownames(hap)))

names(hapInfo) <- c("index","haplotype")
head(hapInfo)
hapInfo$names <- names$IDs[hapInfo$index]
head(hapInfo)
names(hapInfo) <- c("index","haplotype", "IDs")

str(hapInfo)

names_with_pops <- merge(hapInfo, 
                         pops, 
                         by.x = "IDs", 
                         by.y = "Lab.ID..")


# Create the frequency table 
table(names_with_pops$Pop.Structure.Location)
#Atlantic Dry Tortuga      NGOMex      WGOMex 
#6           8           8           2 
pop_freq <- table(names_with_pops$hap, names_with_pops$Pop.Structure.Location)


net <- haploNet(hap)


pdf("../figures/mito_network_males.pdf", h=4, w=4)

plot(net, 
     size = attr(net, "freq"),
     pie = pop_freq,
     scale.ratio=20, cex=0.7,
     show.mutation=1, 
     bg=c("#eac435","#557fc3", "#03cea4", "#fb4d3d"))

legend("bottomleft", colnames(pop_freq), 
       col=c("#eac435","#557fc3", "#03cea4", "#fb4d3d"), pch=19, ncol=1, cex=0.7)

dev.off()

# this plot is very ugly. run it in popart instead and use that for ms. Results are very similar. 



# pop gen ---------------------------------------------------------------

library(adegenet)
library(hierfstat)
dat.gid <- DNAbin2genind(female_seqs)

strata(dat.gid) <- data.frame(Pop.Structure.Location = names_with_pops[,c(4)])
setPop(dat.gid) <- ~Pop.Structure.Location
dat.gid


dat.gid@tab
out <- genind2hierfstat(dat.gid)

pairwise.neifst(out,diploid=FALSE) # gives same as above
a <- genet.dist(out, diploid=FALSE, method = "WC84")
a

# Calculate Nei's GST (FST)
nei.fst <- mmod::pairwise_Gst_Nei(dat.gid)

# command below is wc87 fst
wc.fst <- pairwise.WCfst(out, diploid=F)
boot.ppfst(dat=out,nboot=1000,quant=c(0.025,0.975),diploid=FALSE)

# boot.ppfst runs bootstrapping. but I actually want to do permutation, mix pop assignments
# this is what snpR definitely does and, I think arlequin both do.
# makes sense to test null of no structure.

# Calculate observed FST
observed_fst <- pairwise.WCfst(out, diploid=F)
#               Atlantic        NGOMex        WGOMex Dry Tortuga
#Atlantic             NA 0.68586387 0.76139188 -0.03846154
#NGOMex       0.68586387         NA 0.06582713  0.65445026
#WGOMex       0.76139188 0.06582713         NA  0.72844037
#Dry Tortuga -0.03846154 0.65445026 0.72844037          NA

npop <- nrow(observed_fst)
p_values <- matrix(NA, npop, npop)

nperm <- 1000
perm_dist <- array(NA, dim=c(npop, npop, nperm))

# Run permutations
for(i in 1:nperm) {
  out_perm <- out
  out_perm$pop <- sample(out$pop)
  perm_fst <- pairwise.WCfst(out_perm, diploid=F)
  perm_fst[perm_fst < 0] <- 0
  perm_dist[,,i] <- perm_fst
}

observed_fst[observed_fst < 0] <- 0

# calc pvals
for(i in 1:npop) {
  for(j in 1:npop) {
    p_values[i,j] <- (sum(perm_dist[i,j,] >= observed_fst[i,j])+1)/(nperm+1)
  }
}

p_values


###--------------------------
# plot
###--------------------------
# convert fst to long

melt_fst <- as_tibble(observed_fst, rownames = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")%>%
  mutate(across(c(Var1, Var2), ~str_replace(., "Dry_Tortuga", "Dry Tortuga")))

melt_fst <- subset(melt_fst, !(Var1 != Var2 & is.na(value)))

melt_fst <- melt_fst %>%
  group_by(value) %>%
  filter(!(duplicated(paste0(pmin(Var1, Var2), pmax(Var1, Var2))) & !is.na(value)))

# convert p-values to long:
rownames(p_values) <- row.names(observed_fst)
colnames(p_values) <- colnames(observed_fst)

melt_fst_pval<- as_tibble(p_values, rownames = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")%>%
  mutate(across(c(Var1, Var2), ~str_replace(., "Dry_Tortuga", "Dry Tortuga")))

melt_fst_pval <- melt_fst_pval%>%
  group_by(value) %>%
  filter(!(duplicated(paste0(pmin(Var1, Var2), pmax(Var1, Var2))) & !is.na(value)))


# flip order for pvals, to put on other triangle:
melt_fst_pval2 <- melt_fst_pval[, c("Var2", "Var1", "value")]
colnames(melt_fst_pval2) <- c("Var1", "Var2", "value")

melt_fst$group <- c("fst")
melt_fst_pval2$group <- c("pval")
fst_plot_dat <- rbind(melt_fst, melt_fst_pval2)

fst_plot_dat$group[melt_fst$Var1 == melt_fst$Var2] <- "diagonal"

fst_plot_dat <- fst_plot_dat %>% 
  distinct()


fst_plot_dat$Var1 <- factor(fst_plot_dat$Var1, c("Atlantic", "Dry Tortuga", "NGOMex", "WGOMex"))
fst_plot_dat$Var2 <- factor(fst_plot_dat$Var2, c("Atlantic", "Dry Tortuga", "NGOMex", "WGOMex"))


fst_plot_dat2 <- fst_plot_dat

# drop p-vals
fst_plot_dat3 <- subset(fst_plot_dat2, !group %in% c( "pval"))
fst_plot_dat4 <- subset(fst_plot_dat2, group %in% c( "pval"))
fst_plot_diags <- subset(fst_plot_dat2, group %in% c( "diagonal"))

# need to flip the order of a couple for plots:
fst_plot_dat3[c(2,3,6), c("Var1","Var2")] <- fst_plot_dat3[c(2,3,6), c("Var2","Var1")]
fst_plot_dat4[c(1,2,4), c("Var1","Var2")] <- fst_plot_dat4[c(1,2,4), c("Var2","Var1")]

#fst_plot_dat3[c(7,9), c("Var1","Var2")] <- fst_plot_dat3[c(7,9), c("Var2","Var1")]
#fst_plot_dat4[c(5,6), c("Var1","Var2")] <- fst_plot_dat4[c(5,6), c("Var2","Var1")]

p_fst <- ggplot(fst_plot_dat3, aes(y = Var2, x = Var1, fill = value)) +
  geom_tile(size=0.2, color="black") +
  geom_text(aes(label = round(value, 2)), na.rm = TRUE, size=3) +
  scale_fill_gradient(limits = c(0, 0.77),
                      low = "#FFCCCC", 
                      high = "firebrick3", 
                      na.value = "white") +  
  theme_bw(base_size = 13) +
  labs(x = NULL, y = NULL, fill =  bquote(F[ST])) +
  geom_tile(data = fst_plot_diags, 
            aes(x = Var1, y = Var2),
            fill="grey",
            size=0.2, color="black") +
  geom_tile(data = fst_plot_dat4, 
            aes(x = Var1, y = Var2),
            fill="white",
            size=0.2, color="black") +
  geom_text(data = fst_plot_dat4, 
            aes(label = round(value, 2)),
            size=3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.position = c(0.95, 0.05),
    #legend.justification = c(1, 0),
    legend.background = element_blank(),  # Remove legend background
    legend.margin = margin(0, 0, 0, 0),  # Remove legend margin
    legend.box.margin = margin(0, 0, 0, -10),  # negative value moves legend left
    legend.key.size = unit(0.5, "lines"),  #  reduced legend key size
    legend.title = element_text(size = 9, margin = margin(b = 0)),  # Reduce bottom margin of title
    legend.text = element_blank(),
    #legend.text = element_text(size = 7, margin = margin(l = 1, r = 0)),  # Reduced legend text size
    legend.spacing.y = unit(0.05, "cm"),  # Reduced spacing between legend elements
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p_fst




#### -----------------------------------------------------------
#### -----------------------------------------------------------
#### -----------------------------------------------------------
# nucleotide diversity

mydna <- read.dna("analysis/mitotyping/Pmac_All_Align_complete_CR_truncated_RADSamples_only-DropAmbig.fasta",
                  format="fasta", as.matrix=TRUE)
female_ids <- pops$Lab.ID..[pops$Sex == "M"]
length(female_ids)

# Subset to females only
female_seqs <- mydna[rownames(mydna) %in% female_ids, ]
dim(female_seqs) 
mydna <- female_seqs


seg.sites(mydna)
seg.sites(hap)

nuc.div(mydna)
# [1] 0.001779346

hap.div(mydna)
#0.7681159

# get CI from bootstrapping
n_reps <- 1000
n_seqs <- nrow(mydna)
results <- numeric(n_reps)

for(i in 1:n_reps) {
  # Sample individuals with replacement
  ind_sample <- sample(1:n_seqs, n_seqs, replace=TRUE)
  boot_seqs <- mydna[ind_sample,]
  results[i] <- nuc.div((boot_seqs))
}

ci <- quantile(results, c(0.025, 0.975))
print(paste("Mean pi, boots:", round(mean(results), 5)))
print(paste("95% CI:", round(ci[1], 5), "-", round(ci[2], 5)))
#[1] "95% CI: 0.00138 - 0.00198"

mt_all <- data.frame(
  Population = "All Populations",
  n_samples = 45,
  nuc_div = nuc.div(mydna),
  mean_boot = mean(results), # I know this is wrong
  ci_lower = ci[1],
  ci_upper = ci[2],
  row.names=NULL)


##### --------------------------
# also separate by population:


metadata <- names_with_pops
n_reps=1000

# Split data by population
populations <- unique(metadata$Pop.Structure.Location)

# make empty df to hold results
results_df <- data.frame(
  Population = character(),
  n_samples = numeric(),
  nuc_div = numeric(),
  mean_boot = numeric(),
  ci_lower = numeric(),
  ci_upper = numeric(),
  stringsAsFactors = FALSE
)

for(pop in populations) {
  pop_ids <- metadata$IDs[metadata$Pop.Structure.Location == pop]
  pop_dna <- mydna[pop_ids,]
  
  pop_nuc_div <- nuc.div(pop_dna)
  n_seqs <- length(pop_ids)
  boot_results <- numeric(n_reps)
  
  for(i in 1:n_reps) {
    ind_sample <- sample(1:n_seqs, n_seqs, replace=TRUE)
    boot_seqs <- pop_dna[ind_sample,]
    boot_results[i] <- nuc.div(boot_seqs)
  }
  
  ci <- quantile(boot_results, c(0.025, 0.975))
  
  results_df <- rbind(results_df, data.frame(
    Population = pop,
    n_samples = n_seqs,
    nuc_div = pop_nuc_div,
    mean_boot = mean(boot_results),
    ci_lower = ci[1],
    ci_upper = ci[2],
    row.names=NULL
  ))
}


results_combined <- rbind(results_df, mt_all)
population_order <- c("Atlantic", "Dry Tortuga","All Populations", "NGOMex", "WGOMex")
results_combined$Population <- factor(results_combined$Population, levels = population_order)
legend_order <- c("All Populations", "Atlantic", "Dry Tortuga", "NGOMex", "WGOMex")
results_combined$population_legend <- factor(results_combined$Population, levels = legend_order)

#Population n_samples      nuc_div   mean_boot     ci_lower    ci_upper population_legend
#1          NGOMex         8 0.0019446522 0.001703702 0.0005983545 0.002243829            NGOMex
#2     Dry Tortuga         8 0.0008975318 0.000789454 0.0000000000 0.001196709       Dry Tortuga
#3        Atlantic         6 0.0021640489 0.001807609 0.0006265271 0.002582897          Atlantic
#4          WGOMex         2 0.0041884817 0.002081675 0.0000000000 0.004188482            WGOMex
#5 All Populations        45 0.0017793459 0.001701195 0.0013771910 0.001984217   All Populations

divplot <- ggplot(results_combined, 
                  aes(x = Population, 
                      y = nuc_div, 
                      fill = population_legend, shape=population_legend)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  geom_point(aes(size = ifelse(Population == "All Populations", "All", "Population"))) +
  scale_size_manual(values = c("All" = 8, "Population" = 4), guide = "none") +
  theme_classic(base_size = 12) +
  ylab("Genetic\ndiversity") +
  scale_color_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_fill_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,21,22,23,24))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(
    breaks = c(0, 0.001, 0.002),
    labels = function(x) format(x, scientific = FALSE, big.mark = "")
  )


divplot

ggsave("figures/mt_diversity_males.png", h=5, w=6)
ggsave("figures/mt_diversity_males.pdf", h=5, w=6)

#drop legend:
divplot_noleg <-divplot + 
  theme(legend.position = "none")


# combine plots:

p_fst2 <- p_fst + theme(plot.margin = margin(t = 0.1,  # Top margin
                                             r = 0.5,  # Right margin
                                             b = 0.2,  # Bottom margin
                                             l = 0.5,
                                             unit="cm")) # Left margin
divplot2 <- divplot_noleg + theme(plot.margin = margin(t = 0.05,  # Top margin
                                                       r = 0.5,  # Right margin
                                                       b = 0.1,  # Bottom margin
                                                       l = 0.5,
                                                       unit="cm")) # Left margin
p_out1 <- ggarrange(p_fst2, divplot2, ncol = 1, nrow = 2, labels=c("B", "C"))

ggsave(file="figures/mito_diversity_males.pdf", p_out1, w=3, h=4)

p_fst2 <- p_fst + theme(plot.margin = margin(t = 0.15,  # Top margin
                                             r = 0.1,  # Right margin
                                             b = 0.1,  # Bottom margin
                                             l = 0.1,
                                             unit="cm")) # Left margin
divplot2 <- divplot_noleg + theme(plot.margin = margin(t = 0.15,  # Top margin
                                                       r = 0.1,  # Right margin
                                                       b = 0.0,  # Bottom margin
                                                       l = 0.3,
                                                       unit="cm")) # Left margin
p_out1 <- ggarrange(p_fst2, divplot2, ncol = 2, nrow = 1, labels=c("B", "C"))
ggsave(file="figures/mito_diversity-2_males.pdf", w=5.5, h=2.5)

