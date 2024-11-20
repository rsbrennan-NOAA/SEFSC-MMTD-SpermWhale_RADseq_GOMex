library(tidyverse)
library(scales)
library(ggbeeswarm)

setwd("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis")

dat <- read.table("relatedness.res", header=T)

boxplot(dat$theta, breaks=50)
boxplot(dat$KING)

ggplot(dat, aes(x=KING, y=theta)) +
  geom_point(alpha=0.5, size=2) +
  theme_bw(base_size=14) +
  geom_hline(yintercept=0.065, color = "dodgerblue3") +
  geom_vline(xintercept=0.0442, color = "firebrick3") +
  geom_hline(yintercept=0.08, color = "dodgerblue3", lty=2) +
  geom_vline(xintercept=0.05, color = "firebrick3", lty=2) 
  
ggsave(file="../figures/relatedness_scatter.png", h=5, w=5)


samplelist <- read_delim("indivs.txt",
                         col_names = c("individual"), delim="\t")

# half sibs = 0.125. 
# full sibs 0.25

# see discussion here: https://academic.oup.com/genetics/article/206/1/105/6064207

# check results
# remember indivs are 0 indexed in ngsadmix
half <- dat[which(dat$theta > 0.065 & dat$KING > 0.0442),] 
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
    # M       F       F       F
# half sibs to remove:
  # Pmac120|Pmac125|Pmac131|Pmac101|Pmac130|Pmac052b|Pmac084
# more distant
  #Pmac117|Pmac096

indivrm <- c("Pmac093", "Pmac100", "Pmac125", "Pmac097", "Pmac120", "Pmac131", "Pmac101", "Pmac130", "Pmac52b", "Pmac84", "Pmac117", "Pmac96")
              # M         F           F           F                      M                    M          M         M         M
write.table(file= "C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis/relatedindivs.txt", 
            data.frame(indivrm), row.name=F, col.names = F, quote = F)


##--------------------------------------------------------------------------
# look at relatedness within and between regions.

pops <- read.csv("../SW_Metadata.csv")


# king threshold for cousins. Harder beyond this- https://www.kingrelatedness.com/manual.shtml#WITHIN
# theta thresholds: full sib: 0.25, half: 0.125, 3rd degree- 0.065
dat$id_a <- NA
dat$id_b <- NA
dat$pop_a <- NA
dat$pop_b <- NA
# remove b's

samplelist$individual <- gsub("b", "", samplelist$individual)
for(i in 1:nrow(dat)){
  dat$id_a[i] <- samplelist$individual[(dat$a[i] + 1)]    
  dat$id_b[i] <- samplelist$individual[(dat$b[i] + 1)]
  # add population
  dat$pop_a[i] <- pops$Pop.Structure.Location[which(pops$Lab.ID.. == dat$id_a[i])]
  dat$pop_b[i] <- pops$Pop.Structure.Location[which(pops$Lab.ID.. == dat$id_b[i])]
}

library(dplyr)

dat <- dat %>%
  mutate(comparison = case_when(
    pop_a == "Atlantic" & pop_b == "Atlantic" ~ "Atlantic",
    pop_a == "NGOMex" & pop_b == "NGOMex" ~ "NGOMex",
    pop_a == "Dry Tortuga" & pop_b == "Dry Tortuga" ~ "Dry Tortuga",
    pop_a == "WGOMex" & pop_b == "WGOMex" ~ "WGOMex",
    (pop_a == "Atlantic" & pop_b == "NGOMex") | (pop_a == "NGOMex" & pop_b == "Atlantic") ~ "Atlantic_NGOMex",
    (pop_a == "Atlantic" & pop_b == "Dry Tortuga") | (pop_a == "Dry Tortuga" & pop_b == "Atlantic") ~ "Atlantic_DryTortuga",
    (pop_a == "Atlantic" & pop_b == "WGOMex") | (pop_a == "WGOMex" & pop_b == "Atlantic") ~ "Atlantic_WGOMex",
    (pop_a == "Dry Tortuga" & pop_b == "WGOMex") | (pop_a == "WGOMex" & pop_b == "Dry Tortuga") ~ "DryTortuga_WGOMex",
    (pop_a == "Dry Tortuga" & pop_b == "NGOMex") | (pop_a == "NGOMex" & pop_b == "Dry Tortuga") ~ "DryTortuga_NGOMex",
    (pop_a == "WGOMex" & pop_b == "NGOMex") | (pop_a == "NGOMex" & pop_b == "WGOMex") ~ "WGOMex_NGOMex",
    
    TRUE ~ "other"  # Default case
  ))

dat$comp_color <- dat$comparison

dat$comp_color <- ifelse(dat$pop_a != dat$pop_b, "Between Population", dat$comparison)


p <- ggplot(dat, aes(x=comparison, y=theta, 
                    fill=comp_color,color=comp_color, shape=comp_color)) + 
  geom_quasirandom() +
  geom_violin(color="black", fill=NA) +
  theme_bw(base_size=14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "top",
        legend.title = element_blank())+
  scale_color_manual(values=c("#E69F00", "grey75","#56B4E9", "#009E73", "#CC79A7"))+
  scale_fill_manual(values=c("#E69F00", "grey75","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,21,22,23,24))+    
  xlab(NULL) +
  guides(color = guide_legend(override.aes = list(size = 5)))
p  

ggsave("../figures/population_relatedness.pdf", p, h=5, w=7)
ggsave("../figures/population_relatedness.png", p, h=5, w=7)

out <- dat[,c("id_a", "id_b", "theta", "KING")]

out.order <- out[order(out$theta, decreasing=T),]
head(out.order)

write.csv(out.order, file="../relatedness.csv", quote=F,row.names=F)


# pull out sex of most related:
half$sex_a <- pops$Sex[match(half$id_a, pops$Lab.ID..)]
half$sex_b <- pops$Sex[match(half$id_b, pops$Lab.ID..)]

half[,c("id_a", "id_b", "theta", "KING", "sex_a", "sex_b")]

library(tidyr)

df1 <- unique(data.frame(ID = c(half$id_a,half$id_b), sex = c(half$sex_a,half$sex_b)))
table(df1$sex)

half[half$theta > 0.15,c("id_a", "id_b", "theta", "KING", "sex_a", "sex_b")]

half[half$theta < 0.15 & half$theta > 0.1 ,c("id_a", "id_b", "theta", "KING", "sex_a", "sex_b")]
half[half$theta < 0.1,c("id_a", "id_b", "theta", "KING", "sex_a", "sex_b")]

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Relationship between distance and relatedness:
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

df_dist <- read.csv("../lc_distances_km.csv")
head(df_dist)
colnames(df_dist) <- c("id_a", "id_b", "distance_km")

dat$id_a == df_dist$id_a & dat$id_b == df_dist$id_b | dat$id_a == df_dist$id_b & dat$id_b == df_dist$id_a

for(i in 1:nrow(dat)){
  dat$dist[i] <- df_dist$distance_km[which(df_dist$id_a == dat$id_a[i] & df_dist$id_b == dat$id_b[i])]
}

dat$dist <- as.numeric(dat$dist)
# one pair is 0, change this to 1 for plotting:
dat$dist[dat$dist == 0] <- 1
dat$comp_color <- factor(dat$comp_colo, levels = c("Between Population", "Atlantic", "Dry Tortuga", "NGOMex", "WGOMex"))
custom_sort <- function(x) {
  factor(x, levels = c("Between Population", unique(x[x != "Between Population"])))
}
df_sorted <- dat[order(custom_sort(dat$comp_color)), ]

p <- ggplot(df_sorted, aes(x=(dist), y=theta,fill=comp_color, color=comp_color, shape=comp_color)) +
  geom_point(size=2) +
  scale_color_manual(values=c("grey75", "#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_fill_manual(values= c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,21,22,23,24)) +
  theme_classic(base_size = 14)+
  theme(legend.title=element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4, shape=c(21,21,22,23,24))))+
  theme(legend.position = c(0.8, 0.8)) +
  xlab("Ocean Distance") +
  ylab("Relatedness")

p
ggsave("../figures/population_relatedness_scatter.pdf", p, h=3, w=5)
ggsave("../figures/population_relatedness_scatter.png", p, h=3, w=5)


############################################
# stats:

#just a regular regression
model_out <- lm(theta ~ log10(dist),data=dat)

summary(model_out)
plot(model_out)

# not that good, but its generally in the ballpark of the results below.

# zero inflated beta regression
library(gamlss)
model_zibeta_1 <- gamlss(theta ~ dist, family = BEZI(), data=dat,
                       control = gamlss.control(n.cyc = 200))
summary(model_zibeta_1)

plot(model_zibeta_1)

model_zibeta <- gamlss(theta ~ dist, family = BEINF0, data=dat,
                       control = gamlss.control(n.cyc = 200))
summary(model_zibeta)
plot(model_zibeta)

# both give same results, but BEZI probably makes more sense. 
# BEZI assumes same process generates zeros and non-zero. 
# BEINFO assumes different process generates 0
