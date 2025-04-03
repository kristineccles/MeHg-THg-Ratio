#####################################################################
# Meta Analysis
# By: Kristin Eccles
# Updated: August 18th, 2023
# Written in R Version 4.3.1
#####################################################################

#load libraries
library(mratios)
library(ggplot2)
library(meta)
library(dplyr)
library(viridis)
library(ggplot2)
library(viridis)
library(sjPlot)
library(factoextra)
library(scales)
library(ggpubr)

#load data
df <- read.csv("thg_mehg_metanalysis.csv", header=TRUE, fileEncoding="latin1")
df2 <- read.csv("ecozone_thg.csv", header=TRUE)

#edit data
df$factor <- paste(df$common_name, df$tissue, sep=" ")

#####################################################################
#### EDA ####
# Summary Table
df$entry_count <-1
df$entry_count<- as.numeric(df$entry_count)
species_summary <- as.data.frame(df) %>%
  group_by(df$common_name, df$tissue) %>%
  summarise(mean_thg = mean(thg_mean_ug_g, na.rm=T), 
            sd_thg=sd(thg_mean_ug_g, na.rm=T),
            sum_thg=sum(thg_n, na.rm=T), 
            mean_mehg = mean(mehg_mean_ug_g, na.rm=T), 
            sd_mehg=sd(mehg_mean_ug_g, na.rm=T), 
            sum_mehg=sum(mehg_n, na.rm=T), 
            mean_ratio=mean(mehg_ratio, na.rm=T),
            sd_ratio=sd(mehg_ratio, na.rm=T),
            n=n())
species_summary <-as.data.frame(species_summary)
species_summary <-subset(species_summary, sum_thg>19)
write.csv(species_summary, "species_summary.csv")

#############################################################################3
#### STEP 1: Get the Ratios and SDs ####

# All Data
all_rom <-metacont(data=df, 
                   n.c = thg_n,
                             mean.c = thg_mean_ug_g,
                             sd.c = thg_sd,
                             n.e = mehg_n,
                             mean.e = mehg_mean_ug_g,
                             sd.e =  mehg_sd,
                             sm = "ROM",
                             backtransf = TRUE)
all_rom
# write results to csv
write.csv(all_rom, "rom_results.csv")

# By Ecozones
ecozones_rom <-metacont(data=df2, 
                   n.c = nreplicates,
                   mean.c = xMERCURY,
                   sd.c = sdMERCURY,
                   n.e = nreplicates,
                   mean.e = xMETHYLHG,
                   sd.e =  sdMETHYLHG,
                   sm = "ROM",
                   backtransf = TRUE)
ecozones_rom
# write results to csv
write.csv(ecozones_rom, "ecozones_rom.csv")

#########################################################
#### STEP 2 ####

# load data 
estimates <- read.csv("final_estimated_ratio.csv", fileEncoding="latin1")

# Estimate ratio means by group #
estimates$species_tissue <- paste(estimates$common_name, estimates$tissue)
estimates <- estimates[with(estimates, order(common_name, tissue)),]

#### Tissue Species  ####
# get summary stats tissue and species
species_summary <- as.data.frame(estimates) %>%
  group_by(estimates$common_name, estimates$tissue) %>%
  summarise(count=sum(mehg_n, na.rm=T))
species_summary <-as.data.frame(species_summary)
write.csv(species_summary, "species_tissue_summary_jan19.csv")

# create study, species, tissue, factor
group_means=metamean(mean= mehg_ratio,
                     sd = mehg_ratio_sd,
                     n = mehg_n,
                     data = estimates, 
                     sm = "MRAW",
                     byvar = species_tissue, 
                     backtransf =TRUE)

group_means_results<- cbind(species_summary, group_means$bylevs, group_means$k.all.w, group_means$TE.random.w, 
                            group_means$lower.random.w, group_means$upper.random.w)
colnames(group_means_results) <- c("common_name", "tissue", "n_samples", "group", "n_studies", "proportion", "lower_CI", "upper_CI")
group_means_removed <- subset(group_means_results, n_samples < 20)
length(group_means_removed$common_name)
group_means_results <- subset(group_means_results, n_samples>19)
group_means_results <-na.omit(group_means_results)
write.csv(group_means_results, "species_tissue_meta.csv")

# forest plot 
#group_means_results <-subset(group_means_results, n.studies >1)
# reorder by mean
group_means_results$order <- factor(group_means_results$group, 
                                  levels = group_means_results$group[order(group_means_results$proportion)])

tissue_species <- ggplot(data = group_means_results, aes(x = proportion, y = order)) + 
  geom_point(size=group_means_results$n_samples/100) + 
  scale_color_viridis(alpha = 1, discrete = TRUE, option= "D", name="Tissue")+
  geom_errorbar(aes(xmin = lower_CI , xmax = upper_CI), width = 0.2) +
  theme_bw()+
  facet_wrap(vars(tissue),   scales = "free_y")+
  labs(x="MeHg: THg Mean Ratio (95% CI)", y="Species Tissue", size =12)+
  theme(legend.text=element_text(size=16), text = element_text(size = 20))+
  geom_vline(xintercept=100, linetype="dashed", color = "black")
tissue_species

#Plot figures with dpi=300
cowplot::ggsave2("split_tissue_species.jpg", tissue_species, limitsize = FALSE,
                   width =20, height = 20, dpi = 300)

tissue_species2 <- ggplot(data = group_means_results, aes(x = proportion, y = order, color = tissue)) + 
  geom_point(size=group_means_results$n_samples/100) + 
  geom_errorbar(aes(xmin = lower_CI , xmax = upper_CI), width = 0.2) +
  theme_bw()+
  labs(x="MeHg: THg Mean Ratio (95% CI)", y="Species/Tissue",color="Tissue", size =12)+
  theme(legend.text=element_text(size=16), text = element_text(size = 20))+
  geom_vline(xintercept=100, linetype="dashed", color = "black")
tissue_species2

#Plot figures with dpi=300
cowplot::ggsave2("all_tissue_species.jpg", tissue_species2, 
                 width = 10, height = 20, dpi = 300)

#####################################################################
#### Tissue by Ecozone ####
eco_estimates <- read.csv("hg_metanalysis_estimates_ecozone.csv")
eco_estimates$tissue_subgroup<- paste(eco_estimates$common_nam, eco_estimates$tissue, eco_estimates$ZONE_NAME)
write.csv(eco_estimates, "eco_estimates.csv")
# get summary stats tissue and species
eco_species_summary <- as.data.frame(eco_estimates) %>%
  group_by(tissue_subgroup) %>%
  summarise(count=sum(mehg_n, na.rm=T), mean_thg=mean(thg_mean_u), mean_hg=mean(mehg_mean_))
eco_species_summary <-as.data.frame(eco_species_summary)
write.csv(eco_species_summary, "eco_species_summary.csv")
eco_means=metamean(mean= mehg_ratio,
                     sd = mehg_ratio_sd,
                     n = mehg_n,
                     data = eco_estimates, 
                     sm = "MRAW",
                     byvar = tissue_subgroup, 
                     backtransf =TRUE)

eco_meta_results<- cbind(eco_means$bylevs, eco_means$k.all.w, 
                         eco_means$TE.random.w, eco_means$lower.random.w, eco_means$upper.random.w)
colnames(eco_meta_results) <- c("group", "n_studies", "proportion", "lower_CI", "upper_CI")

#eco_meta_results <- left_join(eco_meta_results,eco_estimates, by=c("grouping1" = "tissue_subgroup"), keep=FALSE)
#eco_meta_results <- subset(eco_meta_results, n_samples >19)
write.csv(eco_meta_results, "eco_meta_results.csv")

#####################################################################
#### Tissue species group ####

estimates$tissue_subgroup<- paste(estimates$grouping, estimates$tissue)
estimates_grouping <- estimates[with(estimates, order(grouping, tissue)),]

# get summary stats tissue 
tissue_summary <- as.data.frame(estimates_grouping) %>%
  group_by(tissue_subgroup) %>%
  summarise(count=sum(mehg_n, na.rm=T))
tissue_summary <-as.data.frame(tissue_summary)

# create study,tissue, species subgroup factor

tissue_means <- metamean(mean= mehg_ratio,
                     sd = mehg_ratio_sd,
                     n = mehg_n,
                     data = estimates_grouping, 
                     sm = "MRAW",
                     byvar = tissue_subgroup, 
                     backtransf =TRUE, 
                     comb.random = TRUE, 
                     comb.fixed = FALSE)

group_tissue_meta<- cbind(tissue_summary, tissue_means$bylevs, tissue_means$k.all.w, 
                          tissue_means$TE.random.w, tissue_means$lower.random.w, tissue_means$upper.random.w)
group_tissue_meta <- subset(group_tissue_meta, count >19)
group_tissue_meta <-na.omit(group_tissue_meta)
colnames(group_tissue_meta) <- c("tissue_group", "n_samples", "group", "n_studies", "proportion", "lower_CI", "upper_CI")
group_means_results <-as.data.frame(group_tissue_meta)
#write.csv(group_tissue_meta, "group_tissue_meta.csv")

# reorder by mean
group_tissue_meta$tissue_order <- factor(group_tissue_meta$tissue_group, 
                                        levels = group_tissue_meta$tissue_group[order(group_tissue_meta$proportion)])
# add tissue group back
unique_group <- estimates[!duplicated(estimates$tissue_subgroup),]
group_tissue_meta<- left_join(group_tissue_meta,  unique_group[,c("grouping","tissue","tissue_subgroup")], 
                              by=c("tissue_group" ="tissue_subgroup"),keep=FALSE)

# subset out fish fat for plotting purposes
group_tissue_meta <- subset(group_tissue_meta, tissue_order != "Fish fat")

#Plot
tissue_species_subgroup <- ggplot(data = group_tissue_meta, aes(x = proportion, y = tissue_order,  color = tissue)) + 
  geom_point(size=group_tissue_meta$n_samples/150) + 
  #scale_color_viridis(alpha = 1, discrete = TRUE, name="Tissue")+
  geom_errorbar(aes(xmin = lower_CI, xmax = upper_CI), width = 0.2) +
  theme_bw()+
  theme(legend.text=element_text(size=16), text = element_text(size = 16))+
  labs(x="MeHg: THg Mean Ratio (95% CI)", y="", color="Tissue", size =12)+
  theme(legend.text=element_text(size=16), text = element_text(size = 20))+
  geom_vline(xintercept=100, linetype="dashed", color = "black")
tissue_species_subgroup
#Plot figures with dpi=300
cowplot::ggsave2("tissue_species_subgroup.jpg", tissue_species_subgroup, 
                 width = 10, height = 10, dpi = 300)



############################################################
# Publication Plot #

combined_plot <- ggarrange(tissue_species2, tissue_species_subgroup,
                               ncol = 2,
                               labels = "AUTO",
                               widths = c(1, 0.75),
                               common.legend = TRUE,
                               legend = "bottom")
combined_plot

# Export as TIFF with 300 DPI
ggsave("combined_plot.tiff", plot = combined_plot, width = 20, height = 15, dpi = 300, units = "in")






