#admixture
library(tidyverse)

setwd("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis")
### what is our most likely k?

# read in CV scores:
cvin <- read.csv("cv.txt", sep=":", header=F)
colnames(cvin) <- c("id", "cv")
# fix the formatting to get K into numeric format
cvin$K <- substr(cvin$id, 4, 4)

# plot the results
p <- ggplot(cvin,aes(x=K,y=cv)) +
  geom_point(size=3)  + geom_line(group=1)

ggsave("../figures/cv.png", p, h=3, w=3)

# actual results:

samplelist <- read_delim("variants_NoLD.fam",
                       col_names = c("individual", "id2", "a", "b", "c", "d"),
                       delim=" ")

# read in all date, in a loop
## first create an empty dataframe
all_data <- tibble(individual=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

# then add all results to this
for (k in 2:8){
  data <- read_delim(paste0("LDthin_numCorrect.",k,".Q"),
                     col_names = paste0("Q",seq(1:k)),
                     delim=" ")
  data$sample <- samplelist$individual
  data$k <- k
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)
}

# add the population label
#all_data$population <- substr(all_data$sample, 1, 2)
#all_data$population <- factor(all_data$population, 
#                              levels=c("GA", "PL", "HP", "BC", "PC", "TR"))

# our orders are off in our vcf. lets re-order these from south to north. 
#orderlist <- read_tsv("~/shared_materials/tutorial_files/population_order.txt",
#                      col_names = "sample")
#all_data$sample<-factor(all_data$sample,levels=orderlist$sample)

all_data %>%
  filter(k == 2) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  #geom_rug(aes(x=sample, y=value, color=population)) +
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",                    
                    labels=c("1","2"))


p <-  
  all_data %>%
  filter(k < 4) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",
                    labels=seq(1:5)) +
  facet_wrap(~k,ncol=1)
p

ggsave("../figures/Admixture_plot.pdf", p, width = 10, height = 9, units="in")




