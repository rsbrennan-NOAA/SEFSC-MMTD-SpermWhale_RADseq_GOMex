

library(vcfR)
library(adegenet)
library(ggpubr)
library(dplyr)

setwd("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis")

vcf <- read.vcfR( "LDthin_numCorrect_nonrelated.vcf", verbose = FALSE )

genl<-vcfR2genlight(vcf)

# add population ids
pops <- read.csv("../SW_Metadata.csv")

genl@ind.names <- gsub("b", "",genl@ind.names)

ids <- data.frame(IDs = genl@ind.names)
result <- ids %>%
  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location), by = c("IDs" = "Lab.ID.."))

genl@pop <- as.factor(result$Pop.Structure.Location)

## Full data set
pca1 <- glPca(genl,center = T, scale = T, nf = 5)
barplot(100*pca1$eig/sum(pca1$eig), col = heat.colors(50), main="PCA Eigenvalues") # retain first 5 axes, incremental decrease after 2
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

#proportion of explained variance by first three axes
a1<-pca1$eig[1]/sum(pca1$eig) # proportion of variation explained by 1st axis
a2<-pca1$eig[2]/sum(pca1$eig) # proportion of variation explained by 2nd axis 
a3<-pca1$eig[3]/sum(pca1$eig) # proportion of variation explained by 3rd axis
pcvar <- data.frame(Axis = c(1:3), Proportion = c(a1,a2,a3))
pcvar

pca1.scores <- as.data.frame(pca1$scores)
pca1.scores$pop <- pop(genl)
pca1.scores$ind <- genl@ind.names

set.seed(89)
num_pops <- length(levels(factor(pca1.scores$pop)))

# plot PC 1 and 2
pca1.p<-ggscatter(pca1.scores, x = "PC1", y = "PC2", color = "pop",
                   ellipse = T, ellipse.level = 0.95, size = 3,
                  xlab = paste0("PC1 (",round(pcvar[1,2]*100,2),"%)"),
                  ylab = paste0("PC2 (",round(pcvar[2,2]*100,2),"%)"))
pca1.p

# #plot  PC 2 and 3: not informative for this dataset
pca2.p<-ggscatter(pca1.scores, x = "PC2", y = "PC3", color = "pop",
                   ellipse = T, ellipse.level = 0.95,
                  xlab = paste0("PC2 (",round(pcvar[2,2]*100,2),"%)"),
                  ylab = paste0("PC3 (",round(pcvar[3,2]*100,2),"%)"))
pca2.p

#ggsave("../Final_plots/Uro_Full_PCA_Axes1_2.svg",pca1.p,width=16,height=16,dpi=600,units="cm")

############ ----------------------------------------------------------------

# follow the k-1 recommendation for the number of PCAs. So k=4
grp <- find.clusters(genl, n.pca=3, max.n.clust=40)

# keep k-1: 3 
# lowest BIC is 11 clusters- plateaus around 11. seems... incorrect.
# do 4 instead.
dapc1 <- dapc(genl, genl@pop, n.pca=3, n.da=4)

#temp <- optim.a.score(dapc1)
col.in <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")
scatter.dapc(dapc1)
scatter.dapc(dapc1, grp=genl@pop)

pdf(file="../figures/dapc.4.pdf", h=4, w=4)
scatter(dapc1, grp=genl@pop, 
        bg="white", pch=c(15,19,17,18), col=col.in,
        cex=2,
        cstar=0,clab=0,
        scree.pca=FALSE,scree.da=FALSE, leg=TRUE, posi.leg="bottomleft")

dev.off()

#### -------------------------------------------------------------------
# try with only 2 groups
grp <- find.clusters(genl, max.n.clust=40)

# keep k-1: 3 
# keep 11 clusters- plateaus around 11.
newpop <- as.character(genl@pop)
newpop <- ifelse(newpop == "Dry Tortuga", "Atlantic",newpop)
newpop <- ifelse(newpop == "NGOMex", "GOMex",newpop)
newpop <- ifelse(newpop == "WGOMex", "GOMex",newpop)
genl2 <- genl
genl2@pop <- as.factor(newpop)
dapc2 <- dapc(genl2, genl2@pop, n.pca=1, n.da=2)

#temp <- optim.a.score(dapc1)

scatter(dapc2)
scatter.dapc(dapc2, grp=genl2@pop)

pdf(file="../figures/dapc.2.pdf", h=4, w=4)
scatter(dapc2, grp=genl2@pop, posi.da="bottomright", 
        bg="white", pch=17:22, col=col.in,
        cex=2,
        cstar=0,clab=0,
        scree.pca=FALSE,scree.da=FALSE, leg=TRUE, posi.leg="topright")
dev.off()



library("poppr")
pramx <- xvalDapc(tab(genl, NA.method = "mean"), pop(genl))

names(pramx)
pramx[-1]

scatter(pramx$DAPC, cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)



###### --------------------------------------------------------------------------
# leave one out, assignment success

# leave out ~30% of samples. 61 indivs. so select 43 to keep (drop 18)
kept.id <- sample(1:nInd(genl), size=43, replace=F)

x <- genl[kept.id,]
x.sup <- genl[!1:nInd(genl) %in% kept.id,]

dapc4 <- dapc(x ,n.pca=3,n.da=4)
pred.sup <- predict.dapc(dapc4, newdata=x.sup )
names(pred.sup)

# assignments of indivs
pred.sup$assign
#coords of scores:
pred.sup$ind.scores
# posterior membership probs
pred.sup$posterior

round(pred.sup$posterior[1:5, 1:4],3)

# plot results:
col <- rainbow(length(levels(pop(x))))
col.points <- transp(col[as.integer(pop(x))],.2)
scatter(dapc4, col=col, bg="white", scree.da=0, pch="",
        cstar=0, clab=0, legend=TRUE)
par(xpd=TRUE)
points(dapc4$ind.coord[,1], dapc4$ind.coord[,2], pch=20,
       col=col.points, cex=5)
col.sup <- col[as.integer(pop(x.sup))]
points(pred.sup$ind.scores[,1], pred.sup$ind.scores[,2], pch=15,
       col=transp(col.sup,.9), cex=2)
add.scatter.eig(dapc4$eig,15,1,2, posi="bottomright", inset=.02)

mean(as.character(pred.sup$assign)==as.character(pop(x.sup)))

table.value(table(pred.sup$assign, pop(x.sup)), col.lab=levels(pop(x.sup)))

# loop to get some estimate of accuracy
n=100
acc.out <- rep(NA, n)

for(i in 1:n){
  # leave out ~30% of samples. 61 indivs. so select 43 to keep (drop 18)
  kept.id <- sample(1:nInd(genl), size=43, replace=F)
  x <- genl[kept.id,]
  x.sup <- genl[!1:nInd(genl) %in% kept.id,]
  # run dapc
  dapc4 <- dapc(x ,n.pca=3,n.da=4)
  pred.sup <- predict.dapc(dapc4, newdata=x.sup )

  # calculate accuracy of assignments. i.e., is anything coming from its actual population
  acc.out[i] <- mean(as.character(pred.sup$assign)==as.character(pop(x.sup)))
  if(i%%50 ==0){print(i)}
}

boxplot(acc.out)
#####
# repeat for k=2

n=100
acc.out.2 <- NA

for(i in 1:n){
  # leave out ~30% of samples. 61 indivs. so select 43 to keep (drop 18)
  kept.id <- sample(1:nInd(genl2), size=43, replace=F)
  x <- genl2[kept.id,]
  x.sup <- genl2[!1:nInd(genl2) %in% kept.id,]
  # run dapc
  dapc4 <- dapc(x ,n.pca=3,n.da=4)
  pred.sup <- predict.dapc(dapc4, newdata=x.sup )
  
  # calculate accuracy of assignments. i.e., is anything coming from its actual population
  acc.out.2[i] <- mean(as.character(pred.sup$assign)==as.character(pop(x.sup)))
  if(i%%50 ==0){print(i)}
}

dplot <- data.frame(k = c(rep("4", length(acc.out)),rep("2", length(acc.out.2))),
          accuracy = c(acc.out,acc.out.2)
            )

library(ggplot2)
library(ggbeeswarm)

p <- ggplot(data=dplot, aes(x=k, y=accuracy)) +
  geom_boxplot(fill="grey60") +
  theme_bw(base_size = 14) +
  ylim(0,1) + 
  xlab("k number")


ggsave(file="../figures/dapc_accuracy.pdf", p, h=5, w=5)
ggsave(file="../figures/dapc_accuracy.png", p, h=5, w=5)

result <- dplot %>%
  group_by(k) %>%
  summarise(
    mean_value = mean(accuracy, na.rm = TRUE),
    sdev = sd(accuracy, na.rm = TRUE),
    se =  sd(accuracy, na.rm = TRUE) / sqrt(n())
  )
result

#  k     mean_value  sdev     se
#  2        0.519   0.134   0.0134
#  4        0.243   0.103   0.0103


