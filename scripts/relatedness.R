library(tidyverse)

setwd("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis")

dat <- read.table("relatedness.res", header=T)

boxplot(dat$theta, breaks=50)
boxplot(dat$KING)

plot(dat$KING, dat$theta)


samplelist <- read_delim("indivs.txt",
                         col_names = c("individual"), delim="\t")

# half sibs = 0.125. 
# full sibs 0.25

# see discussion here: https://academic.oup.com/genetics/article/206/1/105/6064207

# check results
# remember indivs are 0 indexed in ngsadmix
half <- dat[which(dat$theta > 0.065 | dat$KING > 0.0442),] 
  # king threshold for cousins. Harder beyond this- https://www.kingrelatedness.com/manual.shtml#WITHIN
  # theta thresholds: full sib: 0.25, half: 0.125, 3rd degree- 0.065
half$id_a <- NA
half$id_b <- NA
for(i in 1:nrow(half)){
  half$id_a[i] <- samplelist$individual[(half$a[i] + 1)]    
  half$id_b[i] <- samplelist$individual[(half$b[i] + 1)]    
}

half[,c("id_a", "id_b", "theta", "KING")]


## be a bit more lenient, to acct for errors in estimations. 
# check results
# remember indivs are 0 indexed in ngsadmix
half <- dat[which(dat$theta > 0.08 | dat$KING > 0.05),] 
# king threshold for cousins. Harder beyond this- https://www.kingrelatedness.com/manual.shtml#WITHIN
# theta thresholds: full sib: 0.25, half: 0.125, 3rd degree- 0.065
half$id_a <- NA
half$id_b <- NA
for(i in 1:nrow(half)){
  half$id_a[i] <- samplelist$individual[(half$a[i] + 1)]    
  half$id_b[i] <- samplelist$individual[(half$b[i] + 1)]    
}

half[,c("id_a", "id_b", "theta", "KING")]


# full sibs to remove:
  # Pmac093|Pmac100|Pmac125|Pmac097
# half sibs to remove:
  # Pmac120|Pmac125|Pmac131|Pmac101|Pmac130|Pmac052b|Pmac084|
# more distant
  #Pmac117|Pmac096

indivrm <- c("Pmac093", "Pmac100", "Pmac125", "Pmac097", "Pmac120", "Pmac125", "Pmac131", "Pmac101", "Pmac130", "Pmac052b", "Pmac084", "Pmac117", "Pmac096")

write.table(file= "C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis/relatedindivs.txt", 
            data.frame(indivrm), row.name=F, col.names = F, quote = F)

