# figure 3, combined results

#--------------------------------------------------------------------------
# pca

library(ggplot2)
library(stringr)
library(ggbeeswarm)
setwd("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/")

dat <- read.table("analysis/variants_NoLD_unrelated_PCA.eigenvec", header=F)
eigenval <- read.table("analysis/variants_NoLD_unrelated_PCA.eigenval", header=F)

####################
# plot the PCA
####################

# rename our columns, just for ease
colnames(dat) <- c("ID", "ID2", "PC1", "PC2", "PC3", "PC4", colnames(dat)[7:ncol(dat)])

# add a population label:
pops <- read.csv("SW_Metadata.csv")

dat$ID <- gsub("b", "",dat$ID)

result <- dat %>%
  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location, Sex), by = c("ID" = "Lab.ID.."))

eigenval <- read.table("analysis/variants_NoLD_unrelated_PCA.eigenval", header=F)
pve <- data.frame(PC=1:20, pve=round(eigenval$V1/sum(eigenval$V1)*100,1))

result$Pop.Structure.Location[result$Pop.Structure.Location == "WGOMex"] <- "W. Gulf"
result$Pop.Structure.Location[result$Pop.Structure.Location == "NGOMex"] <- "N. Gulf"

# with points
p_pca <- ggplot(result, aes(PC1, PC2, fill=Pop.Structure.Location, shape=Pop.Structure.Location)) +
  geom_point(size =3, color="black") +
  xlab(paste0("PC1: ",pve$pve[1],"% variance")) +
  ylab(paste0("PC2: ",pve$pve[2],"% variance")) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,22,23,24))+
  theme(legend.position = "right",
        legend.title=element_blank()) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white"),
    legend.margin = margin(5, 5, 5, 5),
    legend.text = element_text(size = 12)
  ) +
  guides(
    fill = guide_legend(nrow = 2, ncol = 2, byrow = TRUE, override.aes = list(size = 5)),
    shape = guide_legend(nrow = 2, ncol = 2, byrow = TRUE, override.aes = list(size = 5))
  )

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# fst

library(snpR)
dat <- snpR::read_vcf("analysis/filtered.final.vcf.gz")
pops <- read.csv("SW_Metadata.csv")

# add meta data information:
## population
colnames(dat) <- gsub("b", "",colnames(dat))
ids <- data.frame(IDs = colnames(dat))
result <- ids %>%
  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location, Sex), by = c("IDs" = "Lab.ID.."))
result$Pop.Structure.Location[result$Pop.Structure.Location =="Dry Tortuga"] <- c("Dry_Tortuga")

sample_meta <- data.frame(pop = result$Pop.Structure.Location, sex=result$Sex)
## order the population
#sample_meta$pop <- factor(sample_meta$pop, levels=c("GA", "HP", "BC", "PC", "TR")) 

# assign meta data to dat
sample.meta(dat) <- sample_meta

# calculate fst between the populations
my.dat <- calc_pairwise_fst(dat, facets="pop", method = "WC", boot = 500)
# the bootstrapping here is actually permutation, mixing pop assignments
fst_pvals <- get.snpR.stats(my.dat, facets = "pop", stats = "fst")$fst.matrix$pop$p
#p1 Dry_Tortuga    NGOMex     WGOMex
#<char>       <num>     <num>      <num>
#  1:    Atlantic   0.3712575 0.1317365 0.11776447
#2: Dry_Tortuga          NA 0.1477046 0.10379242
#3:      NGOMex          NA        NA 0.06586826


# convert p-values to long:

dat_fst_pvals2 <- as.matrix(fst_pvals[,2:ncol(fst_pvals)])
rownames(dat_fst_pvals2) <- fst_pvals$p1

melt_fst_pval<- as_tibble(dat_fst_pvals2, rownames = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")%>%
  mutate(across(c(Var1, Var2), ~str_replace(., "Dry_Tortuga", "Dry Tortuga")))

# remove na that isn't diag
melt_fst_pval <- subset(melt_fst_pval, !(Var1 != Var2 & is.na(value)))


# save p-values
write.csv(melt_fst_pval, file="fst_pvalues.csv", quote=F, row.names = F)


# convert fst to long
dat_fst <- get.snpR.stats(my.dat, facets = "pop", stats = "fst")$fst.matrix$pop$fst

dat_fst2 <- as.matrix(dat_fst[,2:ncol(dat_fst)])
rownames(dat_fst2) <- dat_fst$p1

melt_fst <- as_tibble(dat_fst2, rownames = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")%>%
  mutate(across(c(Var1, Var2), ~str_replace(., "Dry_Tortuga", "Dry Tortuga")))

melt_fst <- subset(melt_fst, !(Var1 != Var2 & is.na(value)))

# save fst
write.csv(melt_fst, file="fst.csv", quote=F, row.names = F)



#p_fst <- ggplot(melt_fst, aes(y = Var2, x = Var1, fill = value)) +
#  geom_tile() +
#  geom_text(aes(label = round(value, 3)), na.rm = TRUE) +
#  scale_fill_gradient( low = "#FFCCCC", high = "firebrick3", na.value = "white") +
#  theme_bw(base_size = 14) +
#  labs(x = NULL, y = NULL, fill =  bquote(F[ST])) +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#  theme(
#    axis.text.x = element_text(angle = 45, hjust = 1),
#    legend.position = c(0.95, 0.05),
#    legend.justification = c(1, 0),
#    legend.background = element_blank(),  # Remove legend background
#    legend.margin = margin(0, 0, 0, 0),  # Remove legend margin
#    legend.key.size = unit(0.5, "lines"),  #  reduced legend key size
#    legend.title = element_text(size = 9, margin = margin(b = 0)),  # Reduce bottom margin of title
#    legend.text = element_text(size = 7),  # Reduced legend text size
#    legend.spacing.y = unit(0.05, "cm")  # Reduced spacing between legend elements
#  )


# flip order for pvals, to put on other triangle:
melt_fst_pval2 <- melt_fst_pval[, c("Var2", "Var1", "value")]
colnames(melt_fst_pval2) <- c("Var1", "Var2", "value")

melt_fst$group <- c("fst")
melt_fst_pval2$group <- c("pval")
fst_plot_dat <- rbind(melt_fst, melt_fst_pval2)

fst_plot_dat$group[melt_fst$Var1 == melt_fst$Var2] <- "diagonal"

fst_plot_dat$Var1 <- as.character(fst_plot_dat$Var1 )
fst_plot_dat$Var2 <- as.character(fst_plot_dat$Var2 )
fst_plot_dat$Var1[fst_plot_dat$Var1 == "NGOMex"] <- "N. Gulf"
fst_plot_dat$Var2[fst_plot_dat$Var2 == "NGOMex"] <- "N. Gulf"
fst_plot_dat$Var1[fst_plot_dat$Var1 == "WGOMex"] <- "W. Gulf"
fst_plot_dat$Var2[fst_plot_dat$Var2 == "WGOMex"] <- "W. Gulf"

fst_plot_dat$Var1 <- factor(fst_plot_dat$Var1, c("Atlantic", "Dry Tortuga", "N. Gulf", "W. Gulf"))
fst_plot_dat$Var2 <- factor(fst_plot_dat$Var2, c("Atlantic", "Dry Tortuga", "N. Gulf", "W. Gulf"))


ggplot(fst_plot_dat, aes()) +
  geom_tile(aes(x = Var1, y = Var2, fill = value),
            size = 0.1
  ) 

fst_plot_dat2 <- fst_plot_dat

# add other diags
fst_plot_dat2 <- rbind(fst_plot_dat2, data.frame(
                          Var1 = c("Atlantic", "W. Gulf"),
                          Var2 = c("Atlantic", "W. Gulf"),
                          value = c(NA, NA),
                          group= c("diagonal", "diagonal")
                    )
)

ggplot(fst_plot_dat2, aes()) +
  geom_tile(aes(x = Var1, y = Var2, fill = value),
            size = 0.1
  ) 


# drop p-vals
fst_plot_dat3 <- subset(fst_plot_dat2, !group %in% c( "pval"))
fst_plot_dat4 <- subset(fst_plot_dat2, group %in% c( "pval"))
fst_plot_diags <- subset(fst_plot_dat2, group %in% c( "diagonal"))


p_fst <- ggplot(fst_plot_dat3, aes(y = Var2, x = Var1, fill = value)) +
  geom_tile(size=0.2, color="black") +
  geom_text(aes(label = round(value, 3)), na.rm = TRUE, size=3) +
  scale_fill_gradient(limits = range(fst_plot_dat3$value, na.rm = TRUE),
                      low = "#FFCCCC", 
                      high = "firebrick3", 
                      na.value = "white") +  
  theme_bw(base_size = 13) +
  labs(x = NULL, y = NULL, fill =  bquote(F[ST])) +
  geom_tile(data = fst_plot_diags, 
            aes(x = Var1, y = Var2),
             fill="grey",
            size=0.2, color="black") +
  geom_tile(data = fst_plot_dat4, 
            aes(x = Var1, y = Var2),
            fill="white",
            size=0.2, color="black") +
  geom_text(data = fst_plot_dat4, 
            aes(label = round(value, 2)),
            size=3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.position = c(0.95, 0.05),
    #legend.justification = c(1, 0),
    legend.background = element_blank(),  # Remove legend background
    legend.margin = margin(0, 0, 0, 0),  # Remove legend margin
    legend.box.margin = margin(0, 0, 0, -10),  # negative value moves legend left
    legend.key.size = unit(0.5, "lines"),  #  reduced legend key size
    legend.title = element_text(size = 9, margin = margin(b = 0)),  # Reduce bottom margin of title
    legend.text = element_blank(),
    #legend.text = element_text(size = 7, margin = margin(l = 1, r = 0)),  # Reduced legend text size
    legend.spacing.y = unit(0.05, "cm"),  # Reduced spacing between legend elements
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
  



#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# DAPC accuracy

dat_dapc <- read.csv("analysis/dapc_training_test.csv")

dat_dapc_sub <- dat_dapc %>%
  filter(k %in% c("4pops.3pcs", "2pops.1pcs"))

p_dapc <- ggplot(data=dat_dapc_sub, aes(x=k, y=accuracy)) +
  geom_violin(fill="grey80") +
  geom_boxplot(fill="grey80", width=0.2) +
  theme_bw(base_size = 12)+
  ylim(0,1) + 
  xlab("Number of populations")+ 
  scale_x_discrete(labels = c("2", "4")) +
  ylab("Assignment\naccuracy") +
  geom_point(data=data.frame(k=c("2pops.1pcs", "4pops.3pcs"), accuracy=c(0.5, 0.25)),
             aes(x=k, y=accuracy, color="Null expectation"),
             size=4, shape=18) +
  scale_color_manual(values = c("Null expectation" = "firebrick3")) +
  theme(legend.position = c(0.6, 0.9),
        legend.background = element_rect(fill = NA),
        legend.title = element_blank()) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

  #guides(color = guide_legend(override.aes = list(size=4, shape=18))) +

  



library(ggpubr)

p_fst2 <- p_fst + theme(plot.margin = margin(t = 0.1,  # Top margin
                                      r = 0.5,  # Right margin
                                      b = 0,  # Bottom margin
                                      l = 0.5,
                                      unit="cm")) # Left margin

p_dapc2 <- p_dapc + theme(plot.margin = margin(t = 0.05,  # Top margin
                                               r = 0.5,  # Right margin
                                               b = 0.5,  # Bottom margin
                                               l = 0.5,
                                               unit="cm")) # Left margin

p_out1 <- ggarrange(p_fst2, p_dapc2, ncol = 1, nrow = 2, labels=c("B", "C"))
p_out <- ggarrange(p_pca,p_out1, ncol=2, widths = c(1, 0.75), labels=c("A", "", ""))
ggsave(file="figures/fig3.pdf", w=7, h=4)
ggsave(file="figures/fig3.png", w=7, h=4)


