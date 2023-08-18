#####################################################################
# Forest Plots for Meta Analysis
# By: Kristin Eccles
# Updated: August 18th, 2023
# Written in R Version 4.3.1
#####################################################################

library(ggplot2)
library(viridis)
library(sjPlot)

#### Tissue  Species Group ####

prop_mean_results <- read.csv("species_tissue_results.csv")
summary(prop_mean_results)
#prop_mean_results <-subset(prop_mean_results, n.studies >1)
# reorder by mean
prop_mean_results$order <- factor(prop_mean_results$ Common.Name , 
                                        levels = prop_mean_results$ Common.Name [order(prop_mean_results$mean)])
unique(prop_mean_results$Tissue)
tissue_species <- ggplot(data = prop_mean_results, aes(x = mean, y = order)) + 
  geom_point(size=prop_mean_results$n.samples/100) + 
  scale_color_viridis(alpha = 1, discrete = TRUE, option= "D", name="Tissue")+
  geom_errorbar(aes(xmin = lower_ci , xmax = upper_ci), width = 0.2) +
  theme_bw()+
  facet_wrap(vars(Tissue),   scales = "free_y")+
  labs(x="MeHg: THg Mean Ratio (95% CI)", y="Tissue", size =12)+
  theme(legend.text=element_text(size=16), text = element_text(size = 20))+
  geom_vline(xintercept=100, linetype="dashed", color = "black")
tissue_species

#Plot figures with dpi=300
cowplot::ggsave2("all_tissue_species.jpg", tissue_species, width = 25, height = 25, dpi = 300)

#########################################################################
#### Tissue Species Subgroup ####

meta_tissue_plot <- na.omit(read.csv("tissue_group_results.csv"))

# reorder by mean
meta_tissue_plot$tissue_order <- factor(meta_tissue_plot$grouping, 
                                        levels = meta_tissue_plot$grouping[order(meta_tissue_plot$mean)])

tissue_species_subgroup <- ggplot(data = meta_tissue_plot, aes(x = mean, y = tissue_order,  color = tissue)) + 
  geom_point(size=meta_tissue_plot$n.Samples/150) + 
  scale_color_viridis(alpha = 1, discrete = TRUE, 
                      #labels = c("Brain", "Egg", "Kidney", "Liver", "Muscle", "Plant Tissue", "Skin", "Whole Body"), 
                      name="Tissue")+
  geom_errorbar(aes(xmin = lower_ci, xmax = upper_ci), width = 0.2) +
  theme_bw()+
  theme(legend.text=element_text(size=16), text = element_text(size = 16))+
  labs(x="MeHg: THg Mean Ratio (95% CI)", y="Tissue", size =12)+
  geom_vline(xintercept=100, linetype="dashed", color = "black")
tissue_species_subgroup
#Plot figures with dpi=300
cowplot::ggsave2("tissue_species_subgroup.jpg", tissue_species_subgroup, width = 10, height = 10, dpi = 300)
