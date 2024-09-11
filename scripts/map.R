# make map

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)


# Get the world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# read in locations
setwd("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/")
dat <- read.csv("SW_Metadata.csv")

# read in samples used for genotyping
samplelist <- read_delim("analysis/variants_NoLD.fam",
                         col_names = c("individual", "id2", "a", "b", "c", "d"),
                         delim=" ")

samplelist$indiv2 <- gsub("b", "",samplelist$individual)

sum(dat$Lab.ID.. %in% samplelist$indiv2)
samplelist$indiv2[!(samplelist$indiv2 %in% dat$Lab.ID..)]

dat_in <- dat[dat$Lab.ID.. %in% samplelist$indiv2,]

length(dat_in$Lab.ID..)

usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

pltpt <- data.frame(
  dat_in$Lab.ID..,
  lat = dat_in$Latitude,
  lon = dat_in$Longitude,
  region = dat_in$Pop.Structure.Location
)

# 
p <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = pltpt, aes(x = lon, y = lat, fill=region), 
             shape= 21, color="black", size = 2.5) +
  scale_fill_manual(values=c("#eac435","#557fc3", "#03cea4", "#fb4d3d"))+
  coord_sf() +
  theme_bw() +
  theme(
    #panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    #axis.text = element_blank(),
    #axis.ticks = element_blank()
    legend.position = "top",
  legend.title=element_blank()) +
  xlab("latitude")+
  ylab("longitude") +
  coord_sf(xlim = c(-100, -62), ylim = c(20, 45), expand = FALSE)
p
  
ggsave("figures/map.pdf", p, h=4, w=5)
ggsave("figures/map.png", p, h=4, w=5)


# add admixture props with pi charts
# NOT USING!!!!!
# but keep bc could be useful later. 

#library(scatterpie)
#library(maps)

#qval<- read_delim("analysis/LDthin_numCorrect.2.Q",
#           col_names = c("Q1", "Q2"))

# add qvals to pltpt

#pltpt2 <- cbind(pltpt, qval)
#pltpt2$total <- pltpt2$Q1 + pltpt2$Q2

#usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

#p <- ggplot() +
#  geom_sf(data = world, fill = "grey90", color = "grey70") +
#  geom_sf(data = usa, fill = NA, color = "grey50") +
#  geom_scatterpie(data = pltpt2, 
#                  aes(x = lon, y = lat), 
#                  cols = c("Q1", "Q2"), 
#                  alpha = 0.8, 
#                  color = "black",
#                  linewidth = 0.1) +
#  coord_sf() +
#  theme_minimal() +
#  theme(
#    panel.background = element_rect(fill = "white"),
#    panel.grid = element_blank(),
    #axis.text = element_blank(),
    #axis.ticks = element_blank()
#    legend.position = "none"  # This line removes the legend
#  ) +
#  coord_sf(xlim = c(-100, -62), ylim = c(20, 45), expand = FALSE) +
#  scale_fill_manual(values = c("Q1" = "dodgerblue3", "Q2" = "firebrick2"))
#p

#ggsave("figures/map.admixture.pdf", p, h=6, w=6)
