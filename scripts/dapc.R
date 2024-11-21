

library(vcfR)
library(adegenet)
library(ggpubr)
library(dplyr)
library(poppr)


setwd("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis")

vcf <- read.vcfR( "freebayes_ldthin_unrelated.vcf.gz", verbose = FALSE )

genl<-vcfR2genlight(vcf)

# add population ids
pops <- read.csv("../SW_Metadata.csv")

genl@ind.names <- gsub("b", "",genl@ind.names)

ids <- data.frame(IDs = genl@ind.names)
result <- ids %>%
  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location), by = c("IDs" = "Lab.ID.."))

genl@pop <- as.factor(result$Pop.Structure.Location)

## Full data set
# re-running pca, basically to make sure consistent again with plink, etc. 
pca1 <- glPca(genl,center = T, scale = T, nf = 5)

png("../figures/dapc_eigen.png", h=4, w=4, units="in", res=300)
barplot(100*pca1$eig/sum(pca1$eig), col = heat.colors(50), main="PCA Eigenvalues") # retain first 5 axes, incremental decrease after 2
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
dev.off()

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
                  ylab = paste0("PC2 (",round(pcvar[2,2]*100,2),"%)"),
                  palette =c("#E69F00","#56B4E9", "#009E73", "#CC79A7")
)
ggsave("../figures/dapc_PCA.png",pca1.p, h=5, w=5)

# #plot  PC 2 and 3: not informative for this dataset
pca2.p<-ggscatter(pca1.scores, x = "PC2", y = "PC3", color = "pop",
                   ellipse = T, ellipse.level = 0.95,
                  xlab = paste0("PC2 (",round(pcvar[2,2]*100,2),"%)"),
                  ylab = paste0("PC3 (",round(pcvar[3,2]*100,2),"%)"))
pca2.p



###------------------------------------------------------------------------
###------------------------------------------------------------------------
###------------------------------------------------------------------------
# run dapc
###------------------------------------------------------------------------
###------------------------------------------------------------------------
###------------------------------------------------------------------------

# follow the k-1 recommendation for the number of PCAs. So k=4
grp <- find.clusters(genl, n.pca=3, max.n.clust=40)

# keep k-1: 3 
# lowest BIC is 11 clusters- plateaus around 11. seems... incorrect.
# do 4 instead, based on biology
dapc1 <- dapc(genl, genl@pop, n.pca=3, n.da=4)

#temp <- optim.a.score(dapc1)
col.in <- c("#E69F00","#56B4E9", "#009E73", "#CC79A7")
scatter.dapc(dapc1)
scatter.dapc(dapc1, grp=genl@pop)

png(file="../figures/dapc.4.png", h=4, w=4, units="in", res=300)
scatter(dapc1, grp=genl@pop, 
        bg="white", pch=c(16,15,18,17), col=col.in,
        cex=2,
        cstar=0,clab=0,
        scree.pca=FALSE,scree.da=FALSE, leg=TRUE, posi.leg="bottomleft")

dev.off()

#### -------------------------------------------------------------------
# try with only 2 groups
grp <- find.clusters(genl, max.n.clust=40)

# keep k-1: 3 pcs
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

png(file="../figures/dapc.2.png", h=4, w=4, units="in", res=300)
scatter(dapc2, grp=genl2@pop, posi.da="bottomright", 
        bg="white", pch=c(16,15,18,17), col=col.in,
        cex=2,
        cstar=0,clab=0,
        scree.pca=FALSE,scree.da=FALSE, leg=TRUE, posi.leg="topright")
dev.off()



###### --------------------------------------------------------------------------
###### --------------------------------------------------------------------------
###### --------------------------------------------------------------------------
###### --------------------------------------------------------------------------
# leave one out, assignment success

# leave out ~30% of samples. 61 indivs. so select 43 to keep (drop 18)
kept.id <- sample(1:nInd(genl), size=43, replace=F)

x <- genl[kept.id,]
x.sup <- genl[!(1:nInd(genl)) %in% kept.id,]

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

# repeat to get some estimate of accuracy
n=100
acc.out.4 <- rep(NA, n)
df.out.4 <- as.data.frame(matrix(nrow=n, ncol=4))
colnames(df.out.4) <- c("Atlantic", "DryTortuga", "NGOMex", "WGOMex")

for(i in 1:n){
  # leave out ~30% of samples. 61 indivs. so select 43 to keep (drop 18)
  kept.id <- sample(1:nInd(genl), size=43, replace=F)
  x <- genl[kept.id,]
  x.sup <- genl[!(1:nInd(genl)) %in% kept.id,]
  # run dapc
  dapc4 <- dapc(x , x@pop,n.pca=3,n.da=4)
  pred.sup <- predict.dapc(dapc4, newdata=x.sup )

  # calculate accuracy of assignments. i.e., is anything coming from its actual population
  acc.out.4[i] <- mean(as.character(pred.sup$assign)==as.character(pop(x.sup)))
  
  # get accuracy by population
  #
  Atlantic.tmp <- which(as.character(pop(x.sup)) == "Atlantic")
  if(length(Atlantic.tmp) > 0){
    df.out.4$Atlantic[i] <- mean(as.character(pred.sup$assign[Atlantic.tmp])==as.character(pop(x.sup[Atlantic.tmp])))
  }
  Tortuga.tmp <- which(as.character(pop(x.sup)) == "Dry Tortuga")
  if(length(Tortuga.tmp) > 0){
    df.out.4$DryTortuga[i] <- mean(as.character(pred.sup$assign[Tortuga.tmp])==as.character(pop(x.sup[Tortuga.tmp])))
  }
  NGOMex.tmp <- which(as.character(pop(x.sup)) == "NGOMex")
  if(length(NGOMex.tmp) > 0){
    df.out.4$NGOMex[i] <- mean(as.character(pred.sup$assign[NGOMex.tmp])==as.character(pop(x.sup[NGOMex.tmp])))
  }
  WGOMex.tmp <- which(as.character(pop(x.sup)) == "WGOMex")
  if(length(WGOMex.tmp) > 0){
    df.out.4$WGOMex[i] <- mean(as.character(pred.sup$assign[WGOMex.tmp])==as.character(pop(x.sup[WGOMex.tmp])))
  }
  if(i%%50 ==0){print(i)}
}

boxplot(acc.out.4)

png(filename = "../figures/dapc.4pcs.accuracyPops.png", h=6, w=6, units="in", res=300)
boxplot(df.out.4, ylab="Proportion correct assignment",
        col=c("#eac435","#557fc3", "#03cea4", "#fb4d3d"))
dev.off()

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
  dapc2 <- dapc(x, x@pop, n.pca=1, n.da=2)
  # run prediction
  pred.sup <- predict.dapc(dapc2, newdata=x.sup )
  
  # calculate accuracy of assignments. i.e., is anything coming from its actual population
  acc.out.2[i] <- mean(as.character(pred.sup$assign)==as.character(pop(x.sup)))
  if(i%%50 ==0){print(i)}
}




####-------------------------------------------------------------------------------
####-------------------------------------------------------------------------------
####-------------------------------------------------------------------------------
# find optimal PCs instead of the k-1 approach
####-------------------------------------------------------------------------------
####-------------------------------------------------------------------------------
####-------------------------------------------------------------------------------


# 4 groups

# Calculate optimal number of PCs 
dapcTemp <- dapc(genl, genl@pop, 
                        n.pca=60, n.da = 3)   # n.da = 4(species) - 1 = 3
ascores <- optim.a.score(dapcTemp, smart = FALSE, n.sim = 50)   # Optimal number of PCs: 17

dapc4 <- dapc(genl, genl@pop, 
            n.pca=17, n.da = 3)

scatter(dapc4)

pramx <- xvalDapc(tab(genl, NA.method = "mean"), pop(genl))

names(pramx)
pramx[-1]


scatter(pramx$DAPC, cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topright", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

png(file="../figures/dapc.4.17pcs.png", h=4, w=4, units="in", res=300)
scatter(dapc4, grp=genl@pop, 
        bg="white", pch=c(15,19,17,18), col=col.in,
        cex=2,
        cstar=0,clab=0,
        scree.pca=FALSE,scree.da=FALSE, leg=TRUE, posi.leg="bottomleft")

dev.off()


# get some estimate of accuracy
n=100
acc.out.4.16pcs <- rep(NA, n)

for(i in 1:n){
  # leave out ~30% of samples. 61 indivs. so select 43 to keep (drop 18)
  kept.id <- sample(1:nInd(genl), size=43, replace=F)
  x <- genl[kept.id,]
  x.sup <- genl[!1:nInd(genl) %in% kept.id,]
  # run dapc
  dapc4 <- dapc(x, x@pop,n.pca=17,n.da=3)
  pred.sup <- predict.dapc(dapc4, newdata=x.sup )
  
  # calculate accuracy of assignments. i.e., is anything coming from its actual population
  acc.out.4.16pcs[i] <- mean(as.character(pred.sup$assign)==as.character(pop(x.sup)))
  if(i%%50 ==0){print(i)}
}

boxplot(acc.out.4.16pcs)


# ------------------------------------------------------------------------------
# k = 2 from infer:

grp <- find.clusters(genl2, max.n.clust=40)
# keep all the pcs: 60
# lowest BIC is 2. use this
dapc1 <- dapc(genl2, grp$grp)
# keep 35- 80% rule
# keep 2 discriminant functions
scatter(dapc1)


# Calculate optimal number of PCs 
dapcTemp <- dapc(genl2, genl2@pop, 
                 n.pca=60, n.da = 1)   # n.da = 4(species) - 1 = 3
ascores <- optim.a.score(dapcTemp, smart = FALSE, n.sim = 50)   # Optimal number of PCs: 5
  # but looks terrible

pramx <- xvalDapc(tab(genl2, NA.method = "mean"), pop(genl2)) # also says 5 pcs

names(pramx)
pramx[-1]


scatter(pramx$DAPC, cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

dapc2 <- dapc(genl2, genl2@pop, 
              n.pca=5, n.da = 1)
scatter(dapc2)

png(file="../figures/dapc.2.5pcs.png", h=4, w=4, units="in", res=300)
scatter(dapc2, grp=genl2@pop, 
        bg="white", pch=c(15,19,17,18), col=col.in,
        cex=2,
        cstar=0,clab=0,
        scree.pca=FALSE,scree.da=FALSE, leg=TRUE, posi.leg="bottomleft")

dev.off()


# loop to get some estimate of accuracy
n=100
acc.out.2.5pcs <- rep(NA, n)

for(i in 1:n){
  # leave out ~30% of samples. 61 indivs. so select 43 to keep (drop 18)
  kept.id <- sample(1:nInd(genl2), size=43, replace=F)
  x <- genl2[kept.id,]
  x.sup <- genl2[!1:nInd(genl2) %in% kept.id,]
  # run dapc
  dapc2 <- dapc(x, x@pop,n.pca=5,n.da=1)
  # run prediction
  pred.sup <- predict.dapc(dapc2, newdata=x.sup )
  
  # calculate accuracy of assignments. i.e., is anything coming from its actual population
  acc.out.2.5pcs[i] <- mean(as.character(pred.sup$assign)==as.character(pop(x.sup)))
  if(i%%50 ==0){print(i)}
}

boxplot(acc.out.2.5pcs)



#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------
# plot results
#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------

dplot <- data.frame(k = c(rep("4pops.3pcs", length(acc.out.4)),
                          rep("2pops.1pcs", length(acc.out.2)),
                          rep("4pops.17pcs", length(acc.out.4.16pcs)),
                          rep("2pops.5pcs", length(acc.out.2.5pcs))
                          ),
                    accuracy = c(acc.out.4,acc.out.2,acc.out.4.16pcs,acc.out.2.5pcs)
                    )
dplot$k <- factor(dplot$k, levels=c("2pops.1pcs", "2pops.5pcs", "4pops.3pcs", "4pops.17pcs"))

write.csv(file="dapc_training_test.csv",dplot, row.names=F)

library(ggplot2)
library(ggbeeswarm)

p <- ggplot(data=dplot, aes(x=k, y=accuracy)) +
  geom_violin(fill="grey60") +
  geom_boxplot(fill="grey60", width=0.2) +
  theme_bw(base_size = 14) +
  ylim(0,1) + 
  xlab("k number")
p

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

#k   #PCs        mean_value   sdev      se
#2  1       0.453 0.103  0.0103 
#2  5       0.547 0.127  0.0127 
#4  3        0.231 0.0818 0.00818
#4 17      0.278 0.104  0.0104 


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# make plots that mirror thia,'s increasing PCs and shifting BIC from negative to positive slope
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# df to hold output
outdf <- as.data.frame(matrix(ncol=3, nrow=0))
colnames(outdf) <- c("PCs", "K", "BIC")
for(i in c(3, 10, 20, 40, 60)){

  foo.BIC <- find.clusters(genl, n.pca=i, max.n.clust=30, choose=F)
  tmpout <- data.frame(PCs = i, 
                       K = as.numeric(gsub("K=","",names(foo.BIC$Kstat))),
                       BIC = foo.BIC$Kstat)
  outdf <- rbind(outdf, tmpout)
}

my_labeller <- function(string) {
  labeled_string <- paste(string, "PCs")
  return(labeled_string)
}

outdf$PCs <- as.factor(outdf$PCs)
p <- ggplot(data=outdf, aes(x=K, y=BIC)) +
  geom_line() +
  geom_point() +
  facet_wrap(vars(PCs),scales="free",
             labeller = labeller(PCs = my_labeller)) +
  theme_bw()
p

ggsave(file="../figures/DAPC_PCs_BIC.png",p, w=6, h=4)



