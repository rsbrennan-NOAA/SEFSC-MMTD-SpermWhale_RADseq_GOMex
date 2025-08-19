# diversity via pixy
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(patchwork)


setwd("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis/diversity/")

dat_pops <- read.csv("50kb_pi.txt", header=T, sep="\t")
dat_all <- read.csv("50kb_all_pi.txt", header=T, sep="\t")
dat.na <- dat_pops[!is.na(dat_pops$avg_pi),]
dat_all.na <- dat_all[!is.na(dat_all$avg_pi),]
nrow(dat.na)
# 51748
nrow(dat_all.na)
# 12937

# use only the assembled chromosomes. 

dat_pops.sub <- dat.na[grep("NC_",dat.na$chromosome),]
dat_all.sub <- dat_all.na[grep("NC_",dat_all.na$chromosome),]

nrow(dat_pops.sub)
table(dat_pops.sub$pop)
# 11675 for each pop
table(dat_all.sub$pop)
# 11675 for all


hist(dat_pops.sub$no_sites)

ggplot(dat_pops.sub, aes(x=pop, y=(avg_pi))) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 0.002))
sum(dat_pops.sub$avg_pi == 0)
# 20306

# find average:
# Note that if the user wishes to combine information across windows 
# (e.g. by averaging) after the fact, they should sum the raw counts 
# and recompute the differences/comparisons ratios, 
# and not take an average of the summary statistics themselves.

sum(dat_pops.sub$count_diffs)/sum(dat_pops.sub$count_comparisons)
sum(dat_all.sub$count_diffs)/sum(dat_all.sub$count_comparisons)
#0.00136339

mean(dat_pops.sub$avg_pi, na.rm=T)
# actually pretty similar.

pi <- dat_pops.sub %>%
  group_by(pop) %>%
  summarise(overall_pi = sum(count_diffs) / sum(count_comparisons))
pi
#  pop         overall_pi
# Atlantic       0.00137
# Dry Tortuga    0.00133
# NGOMex         0.00140
# WGOMex         0.00132



#####
# bootstrap CI

# Load required libraries
library(dplyr)
library(boot)
library(purrr)

# function to calculate pi
calculate_pi <- function(data, indices) {
  d <- data[indices,] # need to specify indices here bc this is how boot samples. 
  pi <- sum(d$count_diffs) / sum(d$count_comparisons)
  return(pi)
}

# Function to run bootstrap 
bootstrap_population <- function(pop_data) {
  set.seed(123) 
  boot_results <- boot(data = pop_data, statistic = calculate_pi, R = 10000)
  ci <- boot.ci(boot_results, type = "perc")
  
  return(list(
    population = unique(pop_data$pop),
    estimate = boot_results$t0,
    ci_lower = ci$percent[4],
    ci_upper = ci$percent[5]
  ))
}

# run bootstrap
results_pops <- dat_pops.sub %>%
  group_by(pop) %>%
  group_split() %>% # create of df for each pop
  map_dfr(bootstrap_population) # apply function to each df in list


# Function to run bootstrap no populations
bootstrap_all <- function(pop_data) {
  set.seed(123) 
  boot_results <- boot(data = pop_data, statistic = calculate_pi, R = 10000)
  ci <- boot.ci(boot_results, type = "perc")
  
  return(list(
    estimate = boot_results$t0,
    ci_lower = ci$percent[4],
    ci_upper = ci$percent[5]
  ))
}

# run bootstrap
results_pops <- dat_pops.sub %>%
  group_by(pop) %>%
  group_split() %>% # create of df for each pop
  map_dfr(bootstrap_population) # apply function to each df in list

results_all <- bootstrap_all(dat_all.sub) %>% 
  as_tibble()

results_pops
results_all

#estimate ci_lower ci_upper
#<dbl>    <dbl>    <dbl>
#  1  0.00136  0.00133  0.00140


#population  estimate ci_lower ci_upper
#<chr>          <dbl>    <dbl>    <dbl>
# Atlantic     0.00137  0.00133  0.00141
# Dry Tortuga  0.00133  0.00130  0.00137
# NGOMex       0.00140  0.00136  0.00143
# WGOMex       0.00132  0.00128  0.00136


# plot

# Add population column to results_all and combine with results_pops
results_combined <- bind_rows(
  results_pops,
  results_all %>% mutate(population = "All Populations")
)

results_combined$population <- gsub("NGOMex", "N. Gulf", results_combined$population)
results_combined$population <- gsub("WGOMex", "W. Gulf", results_combined$population)

population_order <- c("Atlantic", "Dry Tortuga","All Populations", "N. Gulf", "W. Gulf")

results_combined$population <- factor(results_combined$population, levels = population_order)
legend_order <- c("All Populations", "Atlantic", "Dry Tortuga", "N. Gulf", "W. Gulf")
results_combined$population_legend <- factor(results_combined$population, levels = legend_order)

results_combined$estimate_10 <- results_combined$estimate * 1000
results_combined$ci_lower_10 <- results_combined$ci_lower * 1000
results_combined$ci_upper_10 <- results_combined$ci_upper * 1000

piplot <- ggplot(results_combined, 
                 aes(x = population, 
                     y = estimate_10, 
                     fill = population_legend, shape=population_legend)) +
  geom_errorbar(aes(ymin = ci_lower_10, ymax = ci_upper_10), width = 0.2) +
  geom_point(aes(size = ifelse(population == "All Populations", "All", "Population"))) +
  scale_size_manual(values = c("All" = 8, "Population" = 4), guide = "none") +
  theme_classic(base_size = 12) +
  ylab(expression(atop("Genetic diversity", paste("(10"^"-3", ")")))) +
  scale_color_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_fill_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,21,22,23,24))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) 


piplot

ggsave("../figures/pi.png", h=5, w=6)





###### -------------------------------------------------------------------------
# re-run with 10kb windows

dat <- read.csv("10kb_pi.txt", header=T, sep="\t")
dat.na <- dat[!is.na(dat$avg_pi),]
nrow(dat.na)
# 63036

# use only the assembled chromosomes. 

dat.sub <- dat.na[grep("NC_",dat.na$chromosome),]

nrow(dat.sub)
table(dat.sub$pop)
# 14107 for each pop

hist(dat.sub$no_sites)

ggplot(dat.sub, aes(x=pop, y=(avg_pi))) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 0.002))
sum(dat.na$avg_pi == 0)
# 30259

# find average:
# Note that if the user wishes to combine information across windows 
# (e.g. by averaging) after the fact, they should sum the raw counts 
# and recompute the differences/comparisons ratios, 
# and not take an average of the summary statistics themselves.

sum(dat.sub$count_diffs)/sum(dat.sub$count_comparisons)
# 0.00136
mean(dat$avg_pi, na.rm=T)
# 0.00138 

pi <- dat.sub %>%
  group_by(pop) %>%
  summarise(overall_pi = sum(count_diffs) / sum(count_comparisons))
pi
#  pop         overall_pi
# Atlantic       0.00137
# Dry Tortuga    0.00133
# NGOMex         0.00140
# WGOMex         0.00132


#####
# bootstrap CI

# Load required libraries
library(dplyr)
library(boot)
library(purrr)

# function to calculate pi
calculate_pi <- function(data, indices) {
  d <- data[indices,] # need to specify indices here bc this is how boot samples. 
  pi <- sum(d$count_diffs) / sum(d$count_comparisons)
  return(pi)
}

# Function to run bootstrap 
bootstrap_population <- function(pop_data) {
  set.seed(123) 
  boot_results <- boot(data = pop_data, statistic = calculate_pi, R = 10000)
  ci <- boot.ci(boot_results, type = "perc")
  
  return(list(
    population = unique(pop_data$pop),
    estimate = boot_results$t0,
    ci_lower = ci$percent[4],
    ci_upper = ci$percent[5]
  ))
}

# run bootstrap
results <- dat.sub %>%
  group_by(pop) %>%
  group_split() %>% # create of df for each pop
  map_dfr(bootstrap_population) # apply function to each df in list


results


#population  estimate ci_lower ci_upper
#<chr>          <dbl>    <dbl>    <dbl>
# Atlantic     0.00137  0.00134  0.00141
# Dry Tortuga  0.00133  0.00130  0.00137
# NGOMex       0.00140  0.00136  0.00143
# WGOMex       0.00132  0.00128  0.00136



#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------


# fis,  obs and exp het, from stacks
setwd("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis/")

in_pops <- read.csv("filtered.final.p.sumstats_summary.tsv", header=T, sep="\t")
in_all <- read.csv("filtered.final.p.sumstats_summary_all.tsv", header=T, sep="\t")
dat <- rbind(in_all, in_pops)
head(dat)
dat$population <- c("All Populations", "Atlantic", "Dry Tortuga", "N. Gulf", "W. Gulf") 
dat$population_legend <- c("All Populations", "Atlantic", "Dry Tortuga", "N. Gulf", "W. Gulf")
population_order <- c("Atlantic", "Dry Tortuga","All Populations", "N. Gulf", "W. Gulf")
dat$population <- factor(dat$population, levels = population_order)


obshet_plot <- ggplot(dat, 
                 aes(x = population, 
                     y = Obs_Het, 
                     fill = population_legend, shape=population_legend)) +
  geom_errorbar(aes(ymin = Obs_Het-StdErr.2, ymax = Obs_Het+StdErr.2), width = 0.2) +
  geom_point(aes(size = ifelse(population == "All Populations", "All", "Population"))) +
  scale_size_manual(values = c("All" = 8, "Population" = 4), guide = "none") +
  theme_classic(base_size = 12) +
  ylab("Observed\nheterozygosity") +
  scale_color_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_fill_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,21,22,23,24))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) 

exphet_plot <- ggplot(dat, 
                      aes(x = population, 
                          y = Exp_Het, 
                          fill = population_legend, shape=population_legend)) +
  geom_errorbar(aes(ymin = Exp_Het-StdErr.3, ymax = Exp_Het+StdErr.3), width = 0.2) +
  geom_point(aes(size = ifelse(population == "All Populations", "All", "Population"))) +
  scale_size_manual(values = c("All" = 8, "Population" = 4), guide = "none") +
  theme_classic(base_size = 12) +
  ylab("Expected\nheterozygosity") +
  scale_color_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_fill_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,21,22,23,24))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) 


# combine obs and expected:

combined_plot <- ggplot(dat, 
                        aes(x = population, 
                            y = Exp_Het, 
                            fill = population_legend, shape=population_legend)) +
  geom_errorbar(aes(ymin = Exp_Het-StdErr.3, ymax = Exp_Het+StdErr.3), width = 0.2) +
  geom_point(aes(size = ifelse(population == "All Populations", "All", "Population"))) +
  geom_errorbar(aes(ymin = Obs_Het-StdErr.2, ymax = Obs_Het+StdErr.2), width = 0.2) +
  geom_point(aes(x = population, 
                 y = Obs_Het,
                 size = ifelse(population == "All Populations", "All", "Population"),
                 color=population_legend),
                 fill="white", stroke=1.5, show.legend = FALSE) +
  scale_size_manual(values = c("All" = 8, "Population" = 4), guide = "none") +
  theme_classic(base_size = 12) +
  ylab("Heterozygosity") +
  scale_color_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_fill_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,21,22,23,24))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  geom_point(data = data.frame(x = c("All Populations", "All Populations"), 
                               y = c(0.255, 
                                     0.25)),
             aes(x = x, y = y), 
             shape = c(21, 21), size = 4, fill = c("black", "white"), color = "black", 
             inherit.aes = FALSE) +
  # Add custom legend text
  geom_text(data = data.frame(x = c("All Populations", "All Populations"), 
                              y = c(0.255, 
                                    0.25),
                              label = c("    Expected", "    Observed")),
            aes(x = x, y = y, label = label), 
            hjust = 0, size = 3, inherit.aes = FALSE) 


combined_plot



#FIS

fis_plot <- ggplot(dat, 
                      aes(x = population, 
                          y = Fis, 
                          fill = population_legend, shape=population_legend)) +
  geom_errorbar(aes(ymin = Fis-StdErr.7, ymax = Fis+StdErr.7), width = 0.2) +
  geom_point(aes(size = ifelse(population == "All Populations", "All", "Population"))) +
  scale_size_manual(values = c("All" = 8, "Population" = 4), guide = "none") +
  theme_classic(base_size = 12) +
  ylab("Inbreeding\ncoefficient") +
  scale_color_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_fill_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,21,22,23,24))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  ylim(0, 0.115)


# Modify each plot to remove x-axis labels except bottom plots
piplot_mod <- piplot + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

obshet_plot_mod <- obshet_plot + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())


# Combine the plots
grid_plots <- (guide_area() / (piplot_mod + obshet_plot_mod) / (exphet_plot + fis_plot)) +
  plot_layout(guides = "collect", heights = (c(0.1, 1,1))) +
  plot_annotation(tag_levels = 'A') &
   theme(legend.position = 'top',plot.tag = element_text(face = 'bold', size=16))

grid_plots

ggsave(file="../figures/fig4.pdf",grid_plots,
       w=8, h=5)
ggsave(file="../figures/fig4.png",grid_plots,
       w=8, h=5)



# 3 panel plot:

ggsave(file="../figures/fig4_revision.pdf",
       ggarrange(piplot, combined_plot, fis_plot, common.legend = T, nrow=1, labels="AUTO"),
       w=7, h=3.3)
ggsave(file="../figures/fig4_revision.png",
       ggarrange(piplot, combined_plot, fis_plot, common.legend = T, nrow=1, labels="AUTO"),
       w=7, h=3.3)




# het pre and post dwh:

dat <- snpR::read_vcf("analysis/filtered.final.vcf.gz")
pops <- read.csv("SW_Metadata.csv")

colnames(dat) <- gsub("b", "",colnames(dat))
ids <- data.frame(IDs = colnames(dat))
result <- ids %>%
  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location, Sex, Collection.Date), by = c("IDs" = "Lab.ID.."))
result$Pop.Structure.Location[result$Pop.Structure.Location =="Dry Tortuga"] <- c("Dry_Tortuga")

sample_meta <- data.frame(pop = result$Pop.Structure.Location, sex=result$Sex,
                          date=result$Collection.Date)

sample_meta[sample_meta$pop == "Atlantic",]
sample_meta[sample_meta$pop == "Dry_Tortuga",]
sample_meta[sample_meta$pop == "NGOMex",]
sample_meta[sample_meta$pop == "WGOMex",]

sample_meta$DWH <- "PreDWH"
sample_meta$date <- as.Date(sample_meta$date, format = "%m/%d/%Y")
sample_meta$DWH[sample_meta$date > as.Date("2003-06-01")] <- "Post-DWH"
table(sample_meta$DWH)
sample_meta <- sample_meta[, -3]
sample.meta(dat) <- sample_meta

# calculate fst between the populations
my.dat <- calc_ho(dat, facets="pop.DWH")
fst_pvals <- get.snpR.stats(my.dat, facets = "pop.DWH", stats = "ho")
fst_pvals

my.dat <- calc_ho(dat, facets="pop")
get.snpR.stats(my.dat, facets = "pop", stats = "ho")

table(paste(sample_meta$pop, sample_meta$DWH))
table(paste( sample_meta$DWH))

# can't do it.


