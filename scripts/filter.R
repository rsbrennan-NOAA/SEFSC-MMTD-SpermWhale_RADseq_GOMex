library(tidyverse)
dat <- read.csv("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/out.imiss", 
                header=T, sep="\t")

hist(dat$F_MISS, breaks=20)

dat[which(dat$F_MISS > 0.4),]

dat <- read.csv("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/filtered.7.depth.ldepth.mean", 
                header=T, sep="\t")

hist(dat$MEAN_DEPTH, breaks=20)

mean(dat$MEAN_DEPTH)
# 17.8
# *3 = 53.4

sum(dat$MEAN_DEPTH > 53.4)
dat[which(dat$MEAN_DEPTH > 53.4),]

# HDplot
library(vcfR)

vcfInput<-read.vcfR("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis/filtered.7_newMAF.vcf.gz")

source("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/scripts/HDplot.R")

HDplotResults<-HDplot(vcfInput)

head(HDplotResults)
mean(HDplotResults$H)

hist(HDplotResults$num_hets)

HDplotResults %>% ggplot()+geom_point(aes(x=H,y=D))

#plot H and ratio
HDplotResults %>% ggplot()+geom_point(aes(x=H,y=ratio))


# sites beyond the distribution are likely paralogs
## D- > than 10 for sure, but maybe even 5?
## H - > than about 0.6


sum((HDplotResults$H > 0.6))
#123
sum((abs(HDplotResults$D) > 6))
#445

# positions to exclude:
datexclude <- HDplotResults[which(HDplotResults$H > 0.6 | abs(HDplotResults$D) > 6),]
posexclude <- datexclude[,1:2]
nrow(posexclude)
#498
write.table(posexclude, file="C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/HD_exclude_newMAF.txt",
            quote=F, col.names = FALSE, row.names=FALSE)

