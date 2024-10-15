library(tidyverse)
library(scales)
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
# half sibs to remove:
  # Pmac120|Pmac125|Pmac131|Pmac101|Pmac130|Pmac052b|Pmac084|
# more distant
  #Pmac117|Pmac096

indivrm <- c("Pmac93", "Pmac100", "Pmac125", "Pmac97", "Pmac120", "Pmac125", "Pmac131", "Pmac101", "Pmac130", "Pmac52b", "Pmac84", "Pmac117", "Pmac96")

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
    
    # Add more conditions as needed
    TRUE ~ "other"  # Default case
  ))

dat$comp_color <- dat$comparison

dat$comp_color <- ifelse(dat$pop_a != dat$pop_b, "between_pop", dat$comparison)

library(ggridges)
ggplot(dat, aes(x=theta, y=comparison, fill=comparison, height = stat(density))) + 
  geom_density_ridges(stat = "binline", bins = 30, scale = 0.95, draw_baseline = FALSE)

library(ggbeeswarm)
p <- ggplot(dat, aes(x=comparison, y=theta, color=comp_color)) + 
  geom_quasirandom() +
  geom_violin(color="black", fill=NA) +
  theme_bw(base_size=14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "top",
        legend.title = element_blank())+
  scale_color_manual(values=c("#eac435","grey75","#557fc3", "#03cea4", "#fb4d3d"))+
    xlab(NULL) 
p  

ggsave("../figures/population_relatedness.pdf", p, h=5, w=7)
ggsave("../figures/population_relatedness.png", p, h=5, w=7)

out <- dat[,c("id_a", "id_b", "theta", "KING")]

out.order <- out[order(out$theta, decreasing=T),]
head(out.order)

write.csv(out.order, file="../relatedness.csv", quote=F,row.names=F)

# ------------------------------------------------------------------------------
# Relationship between distance and relatedness:

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

ggplot(dat, aes(x=(dist), y=(theta))) +
  geom_point() +
  geom_smooth(method="lm")

ggplot(dat, aes(x=log10(dist), y=(theta))) +
  geom_point() +
  geom_smooth(method="lm")

dat[(!is.finite(log10(as.numeric(dat$dist)))),]

# exponential decay
theta.0 <- min(dat$theta) * 0.5  

# Estimate the rest parameters using a linear model
model.0 <- lm((theta) ~ log(dist), data=dat)  
alpha.0 <- exp(coef(model.0)[1])
beta.0 <- coef(model.0)[2]
dat$relatedness <- dat$theta

# Starting parameters
start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
start

model <- nls(relatedness ~ alpha * exp(beta * dist) + theta , data = dat, start = start)
plot(dat$dist, dat$relatedness)
lines(dat$dist, predict(model, list(x = dat$dist)), col = 'skyblue', lwd = 3)
  