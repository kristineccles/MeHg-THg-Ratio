#####################################################################
# Meta Analysis
# By: Kristin Eccles
# Updated: Dec 9th, 2025
# Written in R Version 4.3.1
# Note: Order to Run 02
#####################################################################

#Load Libraries
library(mratios)
library(ggplot2)
library(meta)
library(dplyr)
library(viridis)
library(sjPlot)
library(factoextra)
library(scales)
library(ggpubr)
library(RColorBrewer)

#####################################################################
# Load Data
df  <- read.csv("thg_mehg_metanalysis.csv", header = TRUE, fileEncoding = "latin1")

df$factor <- paste(df$common_name, df$tissue, sep = " ")


tissue_colors <- c(
  "muscle"        = "#4477AA",  # muted blue
  "liver"         = "#44AA99",  # teal-green
  "brain"         = "#AA3377",  # muted magenta
  "kidney"        = "#88CCEE",  # light cyan
  "fat"           = "#DDCC77",  # soft yellow
  "egg"           = "#EEAA33",  # amber
  "heart"         = "#332288",  # deep indigo
  "skin"          = "#999933",  # olive-gray
  "whole body"    = "#777777",  # neutral gray
  "plant"  = "#444444"   # dark gray
)

#####################################################################
#### EDA ####
# Summary Table
df$entry_count  <- 1
df$entry_count  <- as.numeric(df$entry_count)

species_summary <- as.data.frame(df) %>%
  group_by(df$common_name, df$tissue) %>%
  summarise(
    mean_thg   = mean(thg_mean_ug_g,   na.rm = TRUE),
    sd_thg     = sd(thg_mean_ug_g,     na.rm = TRUE),
    sum_thg    = sum(thg_n,            na.rm = TRUE),
    mean_mehg  = mean(mehg_mean_ug_g,  na.rm = TRUE),
    sd_mehg    = sd(mehg_mean_ug_g,    na.rm = TRUE),
    sum_mehg   = sum(mehg_n,           na.rm = TRUE),
    mean_ratio = mean(mehg_ratio,      na.rm = TRUE),
    sd_ratio   = sd(mehg_ratio,        na.rm = TRUE),
    n          = n(),
    wavg_thg = weighted.mean(thg_mean_ug_g, thg_n),
    wavg_mehg = weighted.mean(mehg_mean_ug_g, mehg_n))

species_summary <- as.data.frame(species_summary)
species_summary <- subset(species_summary, sum_thg > 19)

write.csv(species_summary, "species_summary.csv")

#############################################################################
#### STEP 1: Get the Ratios and SDs for each entry ####

# All Data
all_rom <- metacont(
  data       = df,
  n.c        = thg_n,
  mean.c     = thg_mean_ug_g,
  sd.c       = thg_sd,
  n.e        = mehg_n,
  mean.e     = mehg_mean_ug_g,
  sd.e       = mehg_sd,
  sm         = "ROM",
  backtransf = TRUE)
all_rom

# write results to csv
write.csv(all_rom, "rom_results.csv")

#########################################################
#### STEP 2: Create the meta analysis ####

# Load data
estimates <- read.csv("final_estimated_ratio.csv", fileEncoding = "latin1")

# Estimate ratio means by group
estimates$species_tissue <- paste(estimates$common_name, estimates$tissue)
estimates <- estimates[with(estimates, order(common_name, tissue)), ]

# Boxplot 
# EDA looking at variability within group 
df_var <- estimates %>%
  group_by(species_tissue) %>%
  mutate(ratio = mehg_ratio,
    n_samples = sum(thg_n, na.rm = TRUE),
    n_studies = n(),
    factor_lab = paste0(species_tissue,
      "\n(n samples = ", n_samples,
      ", n studies = ", n_studies, ")")) %>%
  filter(max(ratio, na.rm = TRUE) - min(ratio, na.rm = TRUE) > 19,
    mean(ratio, na.rm = TRUE) < 101, 
    n_samples >19,
    n_studies >2)%>%
  ungroup()

ratio_boxplot <- ggplot(
  df_var,
  aes(x = factor_lab, y = ratio, fill = tissue)
) +
  geom_boxplot(
    width = 0.7,
    outlier.alpha = 0.4,
    outlier.size = 1.5,
    position = position_dodge(width = 0.8)
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text  = element_text(size = 15),
    legend.title = element_text(size = 15),
    text         = element_text(size = 15),
    axis.text.y  = element_text(size = 15),
    plot.margin  = margin(t = 10, r = 10, b = 40, l = 60)) +
  labs(
    x = "Species / Tissue",
    y = "MeHg:tHg (%)",
    fill = "Tissue"
  ) +
  scale_fill_manual(
    values = tissue_colors,
    drop = TRUE
  )

ratio_boxplot

# Plot figures with dpi = 300
cowplot::ggsave2(
  filename = "ratio_boxplot.jpg",
  plot     = ratio_boxplot,
  width    = 20,
  height   = 10,
  dpi      = 300
)

#### Tissue & Species ####

# Get summary stats tissue and species
species_summary <- as.data.frame(estimates) %>%
  group_by(estimates$common_name, estimates$tissue) %>%
  summarise(count = sum(mehg_n, na.rm = TRUE))

species_summary <- as.data.frame(species_summary)

write.csv(species_summary, "species_tissue_summary_jan19.csv")

# Create study, species, tissue, factor
group_means <- metamean(
  mean       = mehg_ratio,
  sd         = mehg_ratio_sd,
  n          = mehg_n,
  data       = estimates,
  sm         = "MRAW",
  byvar      = species_tissue,
  backtransf = TRUE)

group_means_results <- cbind(
  species_summary,
  group_means$bylevs,
  group_means$k.all.w,
  group_means$TE.random.w,
  group_means$lower.random.w,
  group_means$upper.random.w)

colnames(group_means_results) <- c(
  "common_name",
  "tissue",
  "n_samples",
  "group",
  "n_studies",
  "proportion",
  "lower_CI",
  "upper_CI")

group_means_removed  <- subset(group_means_results, n_samples < 20)
length(group_means_removed$common_name)

group_means_results  <- subset(group_means_results, n_samples > 19)
group_means_results  <- na.omit(group_means_results)

write.csv(group_means_results, "species_tissue_meta.csv")

# Forest-style Plot (by tissue/species)
group_means_labeled <- group_means_results %>%
  arrange(proportion) %>%
  mutate(
    y_lab = paste0(
      group,
      "(n samples = ", n_samples,
      ", n studies = ", n_studies, ")"
    ),
    y_lab = factor(y_lab, levels = y_lab)
  )

tissue_species <- ggplot(
  data = group_means_labeled,
  aes(x = proportion, y = y_lab, color = tissue)
) +
  geom_point(aes(size = n_samples)) +
  geom_errorbar(
    aes(xmin = lower_CI, xmax = upper_CI),
    width = 0.2
  ) +
  geom_vline(
    xintercept = 100,
    linetype = "dashed",
    color = "black"
  ) +
  theme_bw() +
  labs(
    x     = "MeHg : THg Mean Ratio (95% CI)",
    y     = "",
    color = "Tissue"
  ) +
  scale_color_manual(
    values = tissue_colors,
    limits = names(tissue_colors),
    drop = FALSE) +
  scale_x_continuous(
    breaks = seq(0,
      max(group_means_labeled$upper_CI, na.rm = TRUE),
      by = 25
    )
  ) +
  scale_size(range = c(2, 8), guide = "none") +
  theme(
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 18),
    text         = element_text(size = 20),
    axis.text.y  = element_text(size = 15)
  )

tissue_species


# Plot figures with dpi = 300
cowplot::ggsave2(
  filename = "all_tissue_species.jpg",
  plot     = tissue_species,
  width    = 10,
  height   = 20,
  dpi      = 300
)

#####################################################################
#### Tissue Species Group ####

estimates$tissue_subgroup <- paste(estimates$grouping, estimates$tissue)
estimates_grouping        <- estimates[with(estimates, order(grouping, tissue)), ]

# Get summary stats tissue
tissue_summary <- as.data.frame(estimates_grouping) %>%
  group_by(tissue_subgroup) %>%
  summarise(count = sum(mehg_n, na.rm = TRUE))

tissue_summary <- as.data.frame(tissue_summary)

# Create study, tissue, species subgroup factor
tissue_means <- metamean(
  mean        = mehg_ratio,
  sd          = mehg_ratio_sd,
  n           = mehg_n,
  data        = estimates_grouping,
  sm          = "MRAW",
  byvar       = tissue_subgroup,
  backtransf  = TRUE,
  comb.random = TRUE,
  comb.fixed  = FALSE)

group_tissue_meta <- cbind(
  tissue_summary,
  tissue_means$bylevs,
  tissue_means$k.all.w,
  tissue_means$TE.random.w,
  tissue_means$lower.random.w,
  tissue_means$upper.random.w)

group_tissue_meta <- subset(group_tissue_meta, count > 19)
group_tissue_meta <- na.omit(group_tissue_meta)

colnames(group_tissue_meta) <- c(
  "tissue_group",
  "n_samples",
  "group",
  "n_studies",
  "proportion",
  "lower_CI",
  "upper_CI")

group_means_results <- as.data.frame(group_tissue_meta)
# write.csv(group_tissue_meta, "group_tissue_meta.csv")

# Reorder by mean
group_tissue_meta$tissue_order <- factor(
  group_tissue_meta$tissue_group,
  levels = group_tissue_meta$tissue_group[order(group_tissue_meta$proportion)])

# Add tissue group back
unique_group <- estimates[!duplicated(estimates$tissue_subgroup), ]

group_tissue_meta <- left_join(
  group_tissue_meta,
  unique_group[, c("grouping", "tissue", "tissue_subgroup")],
  by   = c("tissue_group" = "tissue_subgroup"),
  keep = FALSE)

# Subset out fish fat for plotting purposes
group_tissue_meta <- subset(group_tissue_meta, tissue_order != "Fish fat")

# Plot
# Compute a nice max for the x-axis
x_max <- max(group_tissue_meta$upper_CI, na.rm = TRUE)

# Add sample size to group label
group_tissue_meta <- group_tissue_meta %>%
  mutate(
    label = paste0(
      tissue_group,
      "(n samples = ", n_samples,
      ", n studies = ", n_studies, ")"
    )
  )

group_tissue_meta <- group_tissue_meta %>%
  mutate(
    tissue = word(tissue_group, -1)
  )

group_tissue_meta <- group_tissue_meta %>%
  mutate(
    tissue = case_when(
      str_detect(tissue_group, "whole body$") ~ "whole body",
      TRUE ~ word(tissue_group, -1)
    ),
    tissue = factor(tissue, levels = names(tissue_colors))
  )

# Reorder using new label rather than tissue_group alone
tissue_species_subgroup <- ggplot(
  data = group_tissue_meta,
  aes(x = proportion,
    y = reorder(label, proportion),
    color = tissue)) +
  geom_point(aes(size = n_samples)) +
  geom_errorbar(
    aes(xmin = lower_CI, xmax = upper_CI),
    width = 0.2) +
  scale_color_manual(
    values = tissue_colors,
    limits = names(tissue_colors),
    drop = FALSE) +
  scale_size(range = c(2, 8), guide = "none") +
  theme_bw() +
  labs(
    x     = "MeHg:tHg Mean Ratio (95% CI)",
    y     = "",
    color = "Tissue") +
  scale_x_continuous(
    breaks = seq(0, max(group_tissue_meta$upper_CI, na.rm = TRUE), by = 25)) +
  geom_vline(xintercept = 100,
    linetype = "dashed",
    color = "black") +
  theme(legend.text  = element_text(size = 16),
    legend.title = element_text(size = 18),
    text         = element_text(size = 20))
tissue_species_subgroup


# Plot figures with dpi = 300
cowplot::ggsave2(
  filename = "tissue_species_subgroup.jpg",
  plot     = tissue_species_subgroup,
  width    = 10,
  height   = 10,
  dpi      = 300)

####################################################################
#### Publication Plot ####

combined_plot <- ggarrange(
  tissue_species,
  tissue_species_subgroup,
  ncol          = 2,
  labels        = "AUTO",
  #widths        = c(1, 0.75),
  common.legend = FALSE,
  legend        = "bottom",
  font.label    = list(size = 20, face = "bold"))

combined_plot

# Export as TIFF with 300 DPI
ggsave(
  filename = "combined_plot.tiff",
  plot     = combined_plot,
  width    = 25,
  height   = 15,
  dpi      = 300)

#############################################################################
