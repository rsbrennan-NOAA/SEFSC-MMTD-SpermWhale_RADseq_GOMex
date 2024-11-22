# make map

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(ggspatial)


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
  geom_point(data = pltpt, aes(x = lon, y = lat, fill=region, shape=region), 
             color="black", size = 2.5) +
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,22,23,24))+
  coord_sf() +
  theme_bw(base_size=14) +
  theme(
    #panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    #axis.text = element_blank(),
    #axis.ticks = element_blank()
    #legend.position = c(1, 0),
    legend.background = element_rect(fill = NA),  # transparent background
    #legend.box.background = element_rect(color = "black", fill = NA),  # only border, no fill
    legend.margin = margin(0, 0, 0, 0),  # removes space around the legend
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0.1, "cm"),  # reduces space between legend and plot
    legend.text = element_text(margin = margin(l = -2)),  # adjust left margin of text
    #legend.justification = c(1, 0),     
    legend.position="top",
    plot.margin = margin(r = 20, l = 5, b = 10, t = 5),  # Added this line
  legend.title=element_blank()) +
  xlab("longitude") +
  ylab("latitude")+
  coord_sf(xlim = c(-100, -65), ylim = c(22, 43), expand = FALSE)+
  annotation_scale(location = "bl", 
                   width_hint = 0.2, 
                   pad_x = unit(0.45, "in"),  # moves scale bar in from right edge
                   pad_y = unit(0.1, "in"),  # moves scale bar up from bottom edge
                   style = "ticks",
                   unit_category = "metric")
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

# ---------------------------------------------------------------------------------
# pairwise distances:
library(marmap)

#--------------#
#
# Calculate least-cost distances across seas
#
#--------------#

# in part from https://github.com/Tom-Jenkins/seascape_rda_tutorial/blob/master/2.Prepare_spatial_data/2.prepare_spatial_data.R

# Import coordinates of sites
dat <- read.csv("../SW_Metadata.csv")

# read in samples used for genotyping
samplelist <- read_delim("variants_NoLD.fam",
                         col_names = c("individual", "id2", "a", "b", "c", "d"),
                         delim=" ")

samplelist$indiv2 <- gsub("b", "",samplelist$individual)

sum(dat$Lab.ID.. %in% samplelist$indiv2)
samplelist$indiv2[!(samplelist$indiv2 %in% dat$Lab.ID..)]

coords <- dat[dat$Lab.ID.. %in% samplelist$indiv2,]

coords.gps = dplyr::select(coords, Longitude , Latitude)
min(coords.gps$Longitude)

str(coords.gps)

# Get bathymetry data from NOAA using marmap package
bathydata = getNOAA.bathy(lon1 = min(coords.gps$Longitude) -1,
                          lon2 = max(coords.gps$Longitude)+1,
                          lat1 = min(coords.gps$Latitude)-1,
                          lat2 = max(coords.gps$Latitude) +1,
                          resolution = 2)

# Get depth of coordinates
region = dat_in$Pop.Structure.Location

depths = cbind(site = coords$Lab.ID.., region=coords$Pop.Structure.Location,
               get.depth(bathydata, coords.gps, locator = FALSE))
depths


# Plot bathymetry data and coordinates
coords <- coords %>%
  mutate(colors = case_when(
    Pop.Structure.Location == "Atlantic" ~ "#eac435",
    Pop.Structure.Location == "Dry Tortuga" ~ "#557fc3",
    Pop.Structure.Location == "NGOMex" ~ "#03cea4",
    Pop.Structure.Location == "WGOMex" ~ "#fb4d3d",
    TRUE ~ "#999999"  # Default color if none of the above match
  ))

plot(bathydata)
points(coords$Lon, coords$Lat, pch = 21, 
       bg =coords$colors, col = "black", cex = 2)


# Create transition object [long run time]
# Author of marmap recommendations:
# Use a minimum depth of -10 to avoid path crossing land masses
# Use a maximum depth of -200 to limit paths to continental shelf
trans1 = trans.mat(bathydata, min.depth = -0.3, max.depth = NULL)
# save(trans1, file = "transition_object.RData")
# load("transition_object.RData")


# Compute least-cost pat  hs [long run time]
lc_paths = lc.dist(trans1, coords.gps, res = "path")
save(lc_paths, file = "least_cost_paths.RData")
load("least_cost_paths.RData")

# Plot paths on a map
# Visually check that no path overlaps land
plot.bathy(bathydata, image= TRUE, land = TRUE, n = 0,
           bpal = list(c(0, max(bathydata), "grey"),
                       c(min(bathydata), 0, "royalblue")))
lapply(lc_paths, lines, col = "orange", lwd = 2, lty = 1)

# Compute least-cost distances (km) matrix
lc_dist = lc.dist(trans1, coords.gps, res = "dist")

# Convert to matrix, rename columns and rows, and export as csv file
lc_mat = as.matrix(lc_dist)
colnames(lc_mat) = as.vector(coords$Lab.ID..)
rownames(lc_mat) = as.vector(coords$Lab.ID..)
lc_mat

# convert to long format.
df <- as.data.frame(lc_mat)
df$row_names <- rownames(lc_mat)
df_long <- df %>%
  pivot_longer(cols = -row_names, 
               names_to = "column", 
               values_to = "distance_km")
colnames(df_long)<- c("indiv_1", "indiv_2", "distance_km")
df_long <- df_long[which(df_long$indiv_1 != df_long$indiv_2),]
head(df_long)

write.csv(df_long, file="../lc_distances_km.csv", row.names=F)
