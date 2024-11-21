# PCA following dropping the unrelated individuals


# first plot the plink result:
library(ggplot2)
library(patchwork)
setwd("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis")
dat <- read.table("variants_NoLD_unrelated_PCA.eigenvec", header=F)
eigenval <- read.table("variants_NoLD_unrelated_PCA.eigenval", header=F)

# first convert to percentage variance explained
pve <- data.frame(PC=1:20, pve=round(eigenval$V1/sum(eigenval$V1)*100,1))

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot the PC's
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") + 
  ylab("% variance") + 
  theme_classic()

a

ggsave("../figures/pc_var_unrelated.png",
       a, w=4, h=4)


####################
# plot the PCA
####################

# rename our columns, just for ease
colnames(dat) <- c("ID", "ID2", "PC1", "PC2", "PC3", "PC4", colnames(dat)[7:ncol(dat)])

# add a population label:
pops <- read.csv("../SW_Metadata.csv")

dat$ID <- gsub("b", "",dat$ID)

result <- dat %>%
  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location, Sex), by = c("ID" = "Lab.ID.."))



# plot the PCA

d <- ggplot(dat, aes(PC1, PC2, label=ID)) +
  #geom_point(size=4.5, shape=21, color="black", fill="grey48") +
  geom_text(size =3) +
  xlab(paste0("PC1: ",pve$pve[1],"% variance")) +
  ylab(paste0("PC2: ",pve$pve[2],"% variance")) +
  theme_bw() +
  ggtitle("plink: PC1, PC2")
#scale_fill_manual(values=c("#68228B",'#B22222',"#CD6090","#87CEFA", "#1874CD"))
d

ggsave("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/figures/PCA_1_2_plink_ids.png",
       d, w=4, h=4)

# with points
d <- ggplot(result, aes(PC1, PC2, fill=Pop.Structure.Location, shape=Pop.Structure.Location)) +
  geom_point(size =3, color="black") +
  xlab(paste0("PC1: ",pve$pve[1],"% variance")) +
  ylab(paste0("PC2: ",pve$pve[2],"% variance")) +
  theme_classic() +
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,22,23,24))+
  theme(legend.position = "top",
        legend.title=element_blank())



ggsave("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/figures/PCA_1_2_plink.png",
       d, w=4, h=4)



d <- ggplot(dat, aes(PC1, PC3, label=ID)) +
  geom_point(size=4.5, shape=21, color="black") +
  geom_text(size =3) +
  xlab(paste0("PC1: ",pve$pve[1],"% variance")) +
  ylab(paste0("PC3: ",pve$pve[3],"% variance")) +
  theme_bw() +
  ggtitle("plink: PC1, PC3")
 scale_fill_manual(values=c("#68228B",'#B22222',"#CD6090","#87CEFA", "#1874CD"))
d

# save as pdf
ggsave("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/figures/PCA_1_3_plink.png",
       d, w=4, h=4)


####---------------------------------------------------------------------------
# snpRelate pca

library(SNPRelate)
library(SeqArray)
pops <- read.csv("../SW_Metadata.csv")

# save females:
write.table(file="../females_only.txt", pops$Lab.ID..[pops$Sex == "F"], quote=F, row.names=F)

# first run with snprelate
# try with missing data allowed. then with no missing data. 
filename = "freebayes_ldthin_unrelated"
filename.gds = paste0(filename, ".gds")
filename.vcf.gz = paste0(filename, ".vcf.gz")
# 1 . Convert VCF to GDS
SeqArray::seqVCF2GDS(vcf.fn = filename.vcf.gz, out.fn = filename.gds, storage.option="ZIP_RA")

gdsin = SeqArray::seqOpen(filename.gds)
print(paste0("The number of SAMPLES in data: ", length(c(SeqArray::seqGetData(gdsin, "sample.id")))))

print(paste0("The number of SNPs in data: ",  length(c(SeqArray::seqGetData(gdsin, "variant.id")))))

summary(m1 <- SeqArray::seqMissing(gdsin, per.variant=TRUE))
summary(m2 <- SeqArray::seqMissing(gdsin, per.variant=FALSE))
samples <- SeqArray::seqGetData(gdsin, "sample.id")
cbind(samples,m2)[order(-m2),]
hist(m2,breaks=50)

pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, 
                               remove.monosnp = T, maf = 0.05) 

eig = pca.out$eigenval[!is.na(pca.out$eigenval)]
barplot(100*eig/sum(eig), main="PCA Eigenvalues")

# scale 
eig <- 100*eig/sum(eig)

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )
dat$IDs <- gsub("b", "",pca.out$sample.id)

result <- dat %>%
  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location, Sex), by = c("IDs" = "Lab.ID.."))

# with ids
d <- ggplot(result, aes(PC1, PC2, label=IDs)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('snpRelate: PC1, PC2')
d

# with points
d <- ggplot(result, aes(PC1, PC2, fill=Pop.Structure.Location)) +
  geom_point(size =3, color="black", shape=21) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  scale_fill_manual(values=c("#eac435","#557fc3", "#03cea4", "#fb4d3d"))+
  theme(legend.position = "top",
        legend.title=element_blank())
d

ggsave("../figures/PCA_1_2_snpRelate_unrelated.png",
       d, w=4, h=4)


# add sex
d <- ggplot(result, aes(PC1, PC2, fill=Pop.Structure.Location, shape=Sex)) +
  geom_point(size =3, color="black") +
  scale_shape_manual(values=c(21,22)) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  scale_fill_manual(values=c("#eac435","#557fc3", "#03cea4", "#fb4d3d"))+
  theme(legend.position = "top",
        legend.title=element_blank())
d

ggsave("../figures/PCA_1_2_snpRelate_unrelated_sex.png",
       d, w=4, h=4)



#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# Only Females
#--------------------------------------------------------------------------------
pops <- read.csv("../SW_Metadata.csv")
pops$ID <- gsub("b", "",pops$Lab.ID..)

filename = "freebayes_ldthin_unrelated"
filename.gds = paste0(filename, ".gds")
filename.vcf.gz = paste0(filename, ".vcf.gz")
# 1 . Convert VCF to GDS
print(paste0("The number of SAMPLES in data: ", length(c(SeqArray::seqGetData(gdsin, "sample.id")))))

print(paste0("The number of SNPs in data: ",  length(c(SeqArray::seqGetData(gdsin, "variant.id")))))

# subset to only females
gdsin
gds_samples <- SeqArray::seqGetData(gdsin, "sample.id")
gds_samples <- gsub("b", "",gds_samples)

sample_idx <- match(gds_samples, pops$ID)  

# Get the matched sex information
matched_sex <- pops$Sex[sample_idx]
SeqArray::seqSetFilter(gdsin, sample.sel=(matched_sex == "F"))

summary(m1 <- SeqArray::seqMissing(gdsin, per.variant=TRUE))
summary(m2 <- SeqArray::seqMissing(gdsin, per.variant=FALSE))
samples <- SeqArray::seqGetData(gdsin, "sample.id")
cbind(samples,m2)[order(-m2),]
hist(m2,breaks=50)

pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, 
                               remove.monosnp = T, maf = 0.05) 

eig = pca.out$eigenval[!is.na(pca.out$eigenval)]
barplot(100*eig/sum(eig), main="PCA Eigenvalues")

# scale 
eig <- 100*eig/sum(eig)

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )
dat$IDs <- gsub("b", "",pca.out$sample.id)

result <- dat %>%
  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location, Sex), by = c("IDs" = "Lab.ID.."))

eig_df <- data.frame(
  PC = factor(1:20, levels=1:20),
  Percentage = (100*eig/sum(eig))[1:20]
)
p1 <- ggplot(eig_df, aes(x=PC, y=Percentage)) +
  geom_bar(stat="identity",fill="grey50", width=0.7) +
  theme_classic() +
  labs(y="Percentage of variance explained",
       x="Principal Component",
       title="Female only PCA") +
  theme(axis.text.x = element_text(angle=45, hjust=1))
p1

# with points
p2 <- ggplot(result, aes(PC1, PC2, fill=Pop.Structure.Location, shape=Pop.Structure.Location)) +
  geom_point(size =3, color="black") +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,22,23,24))+
  theme(legend.position = "right",
        legend.title=element_blank())
p2


combined_plot <- wrap_plots(p1, p2, widths = c(1, 1))
combined_plot
#ggsave("../figures/PCA_1_2_snpRelate_unrelated.png",
#       d, w=4, h=4)

ggsave("../figures/combined_pca_plot_females.pdf", combined_plot, width=9, height=4)
ggsave("../figures/combined_pca_plot_females.png", combined_plot, width=9, height=4)




#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# drop related indivs, then ld prune. 
# to make sure this doesn't screw anything up.
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------


# plink

dat <- read.table("variants_unrelated_NoLD_2_PCA.eigenvec", header=F)
eigenval <- read.table("variants_unrelated_NoLD_2_PCA.eigenval", header=F)

# first convert to percentage variance explained
pve <- data.frame(PC=1:20, pve=round(eigenval$V1/sum(eigenval$V1)*100,1))

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot the PC's
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") + 
  ylab("% variance") + 
  theme_classic()

a


####################
# plot the PCA
####################

# rename our columns, just for ease
colnames(dat) <- c("ID", "ID2", "PC1", "PC2", "PC3", "PC4", colnames(dat)[7:ncol(dat)])

# add a population label:
pops <- read.csv("../SW_Metadata.csv")

dat$ID <- gsub("b", "",dat$ID)

result <- dat %>%
  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location, Sex), by = c("ID" = "Lab.ID.."))



# plot the PCA

d <- ggplot(dat, aes(PC1, PC2, label=ID)) +
  #geom_point(size=4.5, shape=21, color="black", fill="grey48") +
  geom_text(size =3) +
  xlab(paste0("PC1: ",pve$pve[1],"% variance")) +
  ylab(paste0("PC2: ",pve$pve[2],"% variance")) +
  theme_bw() +
  ggtitle("plink: PC1, PC2")
#scale_fill_manual(values=c("#68228B",'#B22222',"#CD6090","#87CEFA", "#1874CD"))
d

#ggsave("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/figures/PCA_1_2_plink_ids.png",
#       d, w=4, h=4)

# with points
d <- ggplot(result, aes(PC1, PC2, fill=Pop.Structure.Location, shape=Pop.Structure.Location)) +
  geom_point(size =3, color="black") +
  xlab(paste0("PC1: ",pve$pve[1],"% variance")) +
  ylab(paste0("PC2: ",pve$pve[2],"% variance")) +
  theme_classic() +
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,22,23,24))+
  theme(legend.position = "top",
        legend.title=element_blank())



#ggsave("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/figures/PCA_1_2_plink.png",
#       d, w=4, h=4)



d <- ggplot(dat, aes(PC1, PC3, label=ID)) +
  geom_point(size=4.5, shape=21, color="black") +
  geom_text(size =3) +
  xlab(paste0("PC1: ",pve$pve[1],"% variance")) +
  ylab(paste0("PC3: ",pve$pve[3],"% variance")) +
  theme_bw() +
  ggtitle("plink: PC1, PC3")
scale_fill_manual(values=c("#68228B",'#B22222',"#CD6090","#87CEFA", "#1874CD"))
d

# save as pdf
#ggsave("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/figures/PCA_1_3_plink.png",
#       d, w=4, h=4)



# ----------------------------------------------
#snpRelate

filename = "filtered.final_ids"
filename.gds = paste0(filename, ".gds")
filename.vcf.gz = paste0(filename, ".vcf.gz")
# 1 . Convert VCF to GDS
#SeqArray::seqVCF2GDS(vcf.fn = filename.vcf.gz, out.fn = filename.gds, storage.option="ZIP_RA")

gdsin = SeqArray::seqOpen(filename.gds)
print(paste0("The number of SAMPLES in data: ", length(c(SeqArray::seqGetData(gdsin, "sample.id")))))

print(paste0("The number of SNPs in data: ",  length(c(SeqArray::seqGetData(gdsin, "variant.id")))))

summary(m1 <- SeqArray::seqMissing(gdsin, per.variant=TRUE))
summary(m2 <- SeqArray::seqMissing(gdsin, per.variant=FALSE))
samples <- SeqArray::seqGetData(gdsin, "sample.id")
cbind(samples,m2)[order(-m2),]
hist(m2,breaks=50)

indivrm <- c("Pmac093", "Pmac100", "Pmac125", "Pmac097", "Pmac120", "Pmac125", "Pmac131", "Pmac101", "Pmac130", "Pmac052b", "Pmac084", "Pmac117", "Pmac096")

keep = samples[which(!samples %in% indivrm)]
length(keep)
snpset <- SNPRelate::snpgdsLDpruning(gdsin, ld.threshold=0.2, autosome.only = F, start.pos="random", 
                                     num.thread=1, remove.monosnp = T, sample.id = keep)  
snpset.id <- unlist(unname(snpset))

# PCA only on SNPs with a minor allele freq greater than 5%
pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, remove.monosnp = T, maf = 0.05,
                               snp.id=snpset.id,
                               sample.id = keep) # filtering for pruned SNPs

eig = pca.out$eigenval[!is.na(pca.out$eigenval)]
barplot(100*eig/sum(eig), main="PCA Eigenvalues")

# scale 
dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )
dat$IDs <- gsub("b", "",pca.out$sample.id)

result <- dat %>%
  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location), by = c("IDs" = "Lab.ID.."))

# with ids
d <- ggplot(result, aes(PC1, PC2, label=IDs)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('snpRelate: PC1, PC2')
d



























