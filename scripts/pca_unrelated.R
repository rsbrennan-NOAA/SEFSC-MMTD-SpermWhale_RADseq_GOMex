# PCA following dropping the unrelated individuals


# first plot the plink result:
library(ggplot2)
setwd("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis")
dat <- read.table("variants_NoLD_unrelated_PCA.eigenvec", header=F)
eigenval <- read.table("variants_NoLD_unrelated_PCA.eigenval", header=F)

# first convert to percentage variance explained
pve <- data.frame(PC=1:20, pve=round(eigenval$V1/sum(eigenval$V1)*100,1))

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot the PC's
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") + 
  ylab("Percentage variance explained") + 
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
#dat$population <- substr(dat$ID, 1,2)

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

ggsave("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/figures/PCA_1_2_plink.png",
       d, w=4, h=4)


d <- ggplot(dat, aes(PC1, PC3, label=ID)) +
  #geom_point(size=4.5, shape=21, color="black", fill="grey48") +
  geom_text(size =3) +
  xlab(paste0("PC1: ",pve$pve[1],"% variance")) +
  ylab(paste0("PC3: ",pve$pve[3],"% variance")) +
  theme_bw() +
  ggtitle("plink: PC1, PC3"samples)
#scale_fill_manual(values=c("#68228B",'#B22222',"#CD6090","#87CEFA", "#1874CD"))
d

# optional, output your pca to a pdf
ggsave("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/figures/PCA_1_3_plink.png",
       d, w=4, h=4)


####---------------------------------------------------------------------------
# snpRelate pca

library(SNPRelate)
library(SeqArray)
pops <- read.csv("../SW_Metadata.csv")

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
  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location), by = c("IDs" = "Lab.ID.."))

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

