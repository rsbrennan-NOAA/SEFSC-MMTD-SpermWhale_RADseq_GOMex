
# make nexus file for popart:
library(ape)

myseqs <- read.dna("analysis/mitotyping/Pmac_All_Align_complete_CR_truncated_RADSamples_only.fasta",format="fasta",as.matrix=FALSE)
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

setwd("~/projects/spermWhaleRad/")

#data <- read.FASTA("analysis/mitotyping/Pmac_All_Align_completeCR_Jan2020-edit.fasta")
data <- read.FASTA("analysis/mitotyping/Pmac_All_Align_complete_CR_truncated_RADSamples_only.fasta")

checkDNA(data)

length(data)  # How many sequences
hap <- haplotype(data, strict = FALSE, trailingGapsAsN = FALSE)
dist.dna(data, "N", p = TRUE)

#h <- factor(attr(hap, "index"))

#dna_sequences <- as.character(data)
#hap <- haplotype(as.DNAbin(dna_sequences))

# Get haplotype frequencies
#h <- haploFreq(as.DNAbin(dna_sequences))

names <- data.frame(IDs = labels(data))


# h <- factor(attr(hap, "index"))
#pop_freq <- table(h, names_with_pops$Pop.Structure.Location)
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


# Create the frequency table for pie charts
table(names_with_pops$Pop.Structure.Location)
pop_freq <- table(names_with_pops$hap, names_with_pops$Pop.Structure.Location)


net <- haploNet(hap)
#net <- mjn(hap)
#(r <- mjn(data))
#plot(r)

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


#--------------------------------------------------------------------------------
# calc pop gen 
library(adegenet)

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


# Calculate Nei's GST (FST)
nei.fst <- pairwise_Gst_Nei(dat.gid)

# Bootstrap
boot.nei <- bootGst(dat.gid, nboot=999, Gst_statistic="Gst_Nei", quiet=TRUE)

# command below is wc87 fst, I think.
wc.fst <- pairwise.WCfst(out, diploid=F)
boot.ppfst(dat=out,nboot=100,quant=c(0.025,0.975),diploid=FALSE)

# run permutation

perm.test.neifst <- function(data, nperm=999) {
  # Calculate observed FST values
  obs.fst <- pairwise.neifst(data, diploid=FALSE)
  
  # Store number of populations
  n.pops <- nrow(obs.fst)
  pop.names <- rownames(obs.fst)
  
  # Matrix to store permuted values
  perm.values <- array(NA, dim=c(n.pops, n.pops, nperm))
  
  # Perform permutations
  for(i in 1:nperm) {
    # Create copy of data
    perm.data <- data
    # Randomly shuffle population assignments
    perm.data$pop <- sample(data$pop, replace=T)
    # Calculate FST for permuted data
    perm.values[,,i] <- pairwise.neifst(perm.data, diploid=FALSE)
  }
  
  # Calculate p-values
  p.values <- matrix(NA, n.pops, n.pops)
  rownames(p.values) <- pop.names
  colnames(p.values) <- pop.names
  
  # For each pair of populations
  for(i in 1:(n.pops-1)) {
    for(j in (i+1):n.pops) {
      # Calculate how many permuted values are >= observed
      p.values[i,j] <- sum(perm.values[i,j,] >= obs.fst[i,j]) / nperm
      p.values[j,i] <- p.values[i,j]  # matrix is symmetric
    }
  }
  
  # Return results
  return(list(
    observed = obs.fst,
    p.values = p.values
  ))
}

# Run the analysis
set.seed(123)  # for reproducibility
results <- perm.test.neifst(out, nperm=999)

# View results
print("Observed FST values:")
print(results$observed)
print("\nP-values:")
print(results$p.values)









#https://grunwaldlab.github.io/poppr/reference/poppr.amova.html
library(poppr)


amova.result <- poppr.amova(dat.gid, ~Pop.Structure.Location)
amova.result

signif   <- randtest(amova.result, nrepet = 999)
plot(signif)
signif
(neinan <- nei.dist(dat.gid))
unique(dat.gid$strata)
#https://popgen.nescent.org/PopDiffSequenceData.html
library(adegenet)
library(adegenet)
library(hierfstat)
data("nancycats")
a <- genet.dist(nancycats, method = "WC84")
out <- genind2hierfstat(nancycats)
a <- genet.dist(out, method = "WC84")

library(hierfstat)
library(adegenet)
is.genind(dat.gid)


bs.nc<-basic.stats(dat.gid)
dat.gid@pop
hierfstat::genet.dist(dat.gid, diploid=FALSE,method="Nei87")
genet.dist(dat.gid[1:73,],diploid=FALSE,method="Dch")


matFst <- genet.dist(nancycats[1:50, ], method = "Nei87")
nc <- genind2hierfstat(nancycats)

hierfstat::genet.dist(nc,method="Nei87")

dat<- genind2hierfstat(dat.gid, pop=dat.gid@pop)
matFst <- pairwise.neifst(genind2hierfstat(dat.gid, pop=dat.gid@pop))
dat.gid@pop


#dat.pegas <- genind2loci(dat.gid)
#print(dat.pegas, details=TRUE)

library("hierfstat")
obj.hf = genind2hierfstat(dat.gid)
genet.dist(obj.hf, method='Nei87', diploid=F)

library(mmod)

mmod::diff_stats(dat.gid)

#beeData_dist <- dist.multidna(beeData, pool = TRUE)
#amova(beeData_dist ~ populations, data = strata(beeData.gid), nperm = 100)

#mmod::Phi_st_Meirmans(dat.gid) # this function calculates overall PhiST, the Fst analog for DNA sequence data
mmod::pairwise_Gst_Nei(dat.gid, linearized = FALSE) # Calculates pairwise Gst. If linearized = TRUE, it calculates 1/(1- Gst)  
mmod::pairwise_Gst_Hedrick(dat.gid, linearized = FALSE)# Calculates pairwise Gst. If linearized = TRUE, it calculates 1/(1- Gst')  
mmod::pairwise_D(dat.gid, linearized = FALSE, hsht_mean = "harmonic") # Calculates pairwise Gst. If linearized = TRUE, it calculates 1/(1- D)  
bs <- mmod::chao_bootstrap(dat.gid, nreps = 100)
mmod::summarise_bootstrap(bs, statistic=mmod::Gst_Nei)     # for Nei's Gst
plot(data, cex = 0.2)


beeData_dist <- dist.multidna(beeData, pool = TRUE)
amova(beeData_dist ~ populations, data = strata(beeData.gid), nperm = 100)

#### -----------------------------------------------------------
# nucleotide diversity

mydna <- read.dna("mitotyping/Pmac_All_Align_complete CR_Jan2020-edit.fasta", format = "fasta")

ss <- sapply(2:80, function(z){
  length(seg.sites(mydna[1:z, ]))
})
plot(2:80, ss, col = "red", xlab = "No of sequences", ylab = "Segregating sites", las = 1)

seg.sites(mydna)
seg.sites(hap)


nuc.div(mydna)
# [1] 0.001571123
hap.div(mydna)
#0.7388167

# these are the same as here: https://www.sciencedirect.com/science/article/pii/S2352485524001373
# https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14698
# low diversity inflates fst: Measures of Divergence Between Populations and the Effect of Forces that Reduce Variability

# use sapply to loop
nd <- sapply(2:80, function(z){
  nuc.div(mydna[1:z, ])
})

plot(2:80, nd, col = "blue", xlab = "No of sequences", ylab = expression(pi), las = 1)

tajima.test(mydna)

Fst(mydna)

