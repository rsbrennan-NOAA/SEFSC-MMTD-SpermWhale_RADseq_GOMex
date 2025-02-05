
# make nexus file for popart:
library(ape)

myseqs <- read.dna("analysis/mitotyping/Pmac_aligned.fasta",format="fasta",as.matrix=FALSE)
myseqs.names <- names(myseqs)
myseqs.nex.fn <- "analysis/mitotyping/Pmac_All.nex" # output in Nexus format
write.nexus.data(as.character(myseqs),myseqs.nex.fn,interleaved=FALSE,gap="-",
                 missing="?",datablock = FALSE)


#### ----------------------------------
library(pegas)
library(ape)
library(seqinr)
library(tidyverse)
library(dplyr)
library(ggpubr)

data <- read.FASTA("analysis/mitotyping/Pmac_All_Align_complete_CR_truncated_RADSamples_only.fasta")

length(data)  # How many sequences
hap <- haplotype(data, strict = FALSE, trailingGapsAsN = FALSE)

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
#NGOMex        0.3421              
#WGOMex        0.4962 0.0014       
#Dry Tortuga   0.0117 0.5013 0.6537

pairwise.neifst(out,diploid=FALSE) # gives same as above
a <- genet.dist(out, diploid=FALSE, method = "Nei87")
a

# Calculate Nei's GST (FST)
nei.fst <- mmod::pairwise_Gst_Nei(dat.gid)

# command below is wc87 fst, I believe
wc.fst <- pairwise.WCfst(out, diploid=F)
boot.ppfst(dat=out,nboot=1000,quant=c(0.025,0.975),diploid=FALSE)

# boot.ppfst runs bootstrapping. but I actually want to do permutation, mix pop assignments
# this is what snpR definitely does and, I think arlequin both do.
# makes sense to test null of no structure.

# Calculate observed FST
observed_fst <- pairwise.WCfst(out, diploid=F)
#               Atlantic        NGOMex        WGOMex Dry Tortuga
#Atlantic            NA  0.3547399965  0.5138892625  0.01360011
#NGOMex      0.35474000            NA -0.0005847247  0.50172523
#WGOMex      0.51388926 -0.0005847247            NA  0.65369306
#Dry Tortuga 0.01360011  0.5017252270  0.6536930562          NA

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
  scale_fill_gradient(limits = range(fst_plot_dat3$value, na.rm = TRUE),
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


### ----------------------
# amova:
library(poppr)

dat.gc <- as.genclone(dat.gid)
amova.result <- poppr.amova(dat.gc, ~Pop.Structure.Location)
amova.test <- randtest(amova.result) # Test for significance
plot(amova.test)

amova.test

set.seed(9001)
dat.gc.new <- dat.gc
strata(dat.gc.new) <- data.frame(Pop = strata(dat.gc)[sample(nInd(dat.gc)),])
dat.gc.new.amova<- poppr.amova(dat.gc.new, ~Pop)
dat.gc.new.amova
dat.new.amova.test<- randtest(dat.gc.new.amova, nrepet = 999)


#### -----------------------------------------------------------
#### -----------------------------------------------------------
#### -----------------------------------------------------------
# nucleotide diversity
#### -----------------------------------------------------------
#### -----------------------------------------------------------
#### -----------------------------------------------------------

mydna <- read.dna("analysis/mitotyping/Pmac_aligned.fasta",format="fasta")

ss <- sapply(2:73, function(z){
  length(seg.sites(mydna[1:z, ]))
})
plot(2:73, ss, col = "red", xlab = "No of sequences", ylab = "Segregating sites", las = 1)

seg.sites(mydna)
seg.sites(hap)


nuc.div(mydna)
# [1] 0.001671559

hap.div(mydna)
#0.7366819

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


mydna <- read.dna("mitotyping/Pmac_All_Align_complete_CR_truncated_RADSamples_only.fasta",format="fasta")



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

ggsave("../figures/mt_diversity.png", h=5, w=6)
ggsave("../figures/mt_diversity.pdf", h=5, w=6)

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

ggsave(file="../figures/mito_diversity.pdf", w=3, h=4)

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
ggsave(file="../figures/mito_diversity-2.pdf", w=5.5, h=2.5)






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




