#####################################################################
# Effect of co-variates on MeHg:THg ratios - Age, Sex, Fork Length
# By: Kristin Eccles
# Updated: Dec 18th, 2025
# Written in R Version 4.3.1
# Note: Can be run independently 
#####################################################################

#load libraries
library(ggplot2)
library(ggpubr)
library(cowplot)
library(sjPlot)
library(corrplot)
library(stats)
library(lmtest)
library(car)
library(multcomp)
library(meta)
library(metafor)
library(ggpmisc)

####################################################################
#### Load Data ####

yp_df <- read.csv("batchelar_yellow_perch.csv")
sex_color <- c("Male" = "#4D79F6", "Female" = "#E66465")

####################################################################
#### Correlation plots (Figure S2): MeHg ratio vs age / weight / length ####

vars <- c("mehg_ratio", "mean_age", "mean_weight_g", "mean_fork_length_mm")

# Clean labels: remove "_" and Title Case
pretty_labels <- vars |>
  str_replace_all("_", " ") |>
  str_to_title()

# Split by sex
yp_male   <- yp_df %>% filter(sex == "Male")
yp_female <- yp_df %>% filter(sex == "Female")

# Site/sample counts
n_site_male    <- nrow(yp_male)
n_sample_male  <- sum(yp_male$thg_n)
n_site_female  <- nrow(yp_female)
n_sample_female <- sum(yp_female$thg_n)

# Correlation matrices
cor_male   <- cor(yp_male[, vars],   use = "pairwise.complete.obs")
cor_female <- cor(yp_female[, vars], use = "pairwise.complete.obs")

colnames(cor_male)   <- rownames(cor_male)   <- pretty_labels
colnames(cor_female) <- rownames(cor_female) <- pretty_labels

# ggplot-style correlation plots
p_male <- ggcorrplot(
  cor_male,
  lab      = TRUE,
  lab_col  = "black",
  colors   = c("red", "white", "blue"),
  ggtheme  = theme_minimal()
) +

  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

p_female <- ggcorrplot(
  cor_female,
  lab      = TRUE,
  lab_col  = "black",
  colors   = c("red", "white", "blue"),
  ggtheme  = theme_minimal()) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

fig_s2 <- ggarrange(
  p_male, p_female,
  ncol          = 2,
  common.legend = TRUE,
  legend        = "bottom",
  labels = "AUTO")

fig_s2

# Export correlation figure
ggsave(
  "Figure_S2_YellowPerch_Male_Female.jpg",
  fig_s2, width = 10, height = 5, dpi = 600)

####################################################################
#### Sex effect: boxplot + + meta-analysis ####

# Boxplot of MeHg ratio by sex 
yp_boxplot <- ggplot(yp_df, aes(x = sex, y = mehg_ratio, fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.3) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.6) +
  scale_fill_manual(values = c("#4D79F6", "#E66465")) +
  labs(x = "Sex", y = "MeHg Ratio (%)", fill = "Sex") +
  theme_bw(base_size = 12)

yp_boxplot

yp_boxplot_ttest <- yp_boxplot +
  stat_compare_means(
    method  = "t.test",
    formula = mehg_ratio ~ sex,
    var.equal = FALSE,          # Welch t-test (unbalanced)
    label  = "p.format",
    label.y = max(yp_df$mehg_ratio) * 1.05)

yp_boxplot_ttest

# Meta-analysis: sex effect on MeHg ratio
m.gen <- metagen(
  TE       = mehg_ratio,
  seTE     = mehg_ratio_se,
  studlab  = first_author,
  data     = yp_df,
  sm       = "SMD",
  fixed    = FALSE,
  random   = TRUE,
  method.tau = "REML")
m.gen

# Subgroup meta-analysis by sex
sex_effect <- update(m.gen, subgroup = sex, tau.common = TRUE)
sex_effect   # no effect of sex in the meta-analysis

# Simple t-test for comparison
t.test(mehg_ratio ~ sex, data = yp_df)  # consistent with meta-analysis

####################################################################
#### Covariate effects: Age, Weight, Fork length ####


## Age
age_plot <- age_plot <- ggplot(yp_df,aes(mean_age, mehg_ratio, color = sex, size = thg_n)) +
  geom_point() +
  scale_color_manual(values = sex_color) +
  scale_size_continuous(
    range = c(2, 10),
    name = "n Samples",
    guide = guide_legend(
      override.aes = list(
        linetype = "blank",   # remove smoother line
        shape = 16,           # round point
        color = "black" ,
        fill = NA )))+
  geom_smooth(method = lm, se = TRUE) +
  stat_poly_eq(
    formula = y~x,
    aes(grp.label = sex,
        label=paste(..grp.label..,..eq.label..,..rr.label..,..p.value.label..,sep="~~~")),
    parse=TRUE, size=5, label.x.npc="right") +

  labs(x="Mean Age (years)", y="MeHg Ratio (%)", color="Sex") +
  theme_bw()
age_plot

age_lm  <- lm(mehg_ratio ~ mean_age + sex, data = yp_df)
summary(age_lm)
Anova(age_lm, type = 3)

age_lm2 <- lm(mehg_ratio ~ mean_age, data = yp_df)
summary(age_lm2)

ancova_age <- aov(mehg_ratio ~ mean_age + sex, data = yp_df)
summary(ancova_age)

## Weight
weight_plot <- ggplot(yp_df,
                      aes(x = mean_weight_g,
                          y = mehg_ratio,
                          color = sex, size = thg_n)) +
  geom_point() +
  scale_color_manual(values = sex_color) +
  scale_size_continuous(
    range = c(2, 10),
    name = "n Samples",
    guide = guide_legend(
      override.aes = list(
        linetype = "blank",   # remove smoother line
        shape = 16,           # round point
        color = "black" ,
        fill = NA )))+
  geom_smooth(method = lm, se = TRUE) +
  stat_poly_eq(
    formula = y~x,
    aes(grp.label = sex,
        label=paste(..grp.label..,..eq.label..,..rr.label..,..p.value.label..,sep="~~~")),
    parse=TRUE, size=5, label.x.npc="right") +
  labs(x="Mean Age (years)", y="MeHg Ratio (%)", color="Sex") +
  theme_bw()
weight_plot

weight_lm  <- lm(mehg_ratio ~ mean_weight_g + sex, data = yp_df)
summary(weight_lm)
Anova(weight_lm, type = 3)

weight_lm2 <- lm(mehg_ratio ~ mean_weight_g, data = yp_df)
summary(weight_lm2)   # not significant

## Fork length
fork_plot <- ggplot(yp_df,
                    aes(x = mean_fork_length_mm,
                        y = mehg_ratio,
                        color = sex, size = thg_n)) +
  geom_point() +
  scale_color_manual(values = sex_color) +
  scale_size_continuous(
    range = c(2, 10),
    name = "n Samples",
    guide = guide_legend(
      override.aes = list(
        linetype = "blank",   # remove smoother line
        shape = 16,           # round point
        color = "black" ,
        fill = NA )))+
  geom_smooth(method = lm, se = TRUE) +
  stat_poly_eq(
    formula = y~x,
    aes(grp.label = sex,
        label=paste(..grp.label..,..eq.label..,..rr.label..,..p.value.label..,sep="~~~")),
    parse=TRUE, size=5, label.x.npc="right") +
  labs(x="Mean Age (years)", y="MeHg Ratio (%)", color="Sex") +
  theme_bw()
fork_plot

fork_lm  <- lm(mehg_ratio ~ mean_fork_length_mm + sex, data = yp_df)
summary(fork_lm)

fork_lm2 <- lm(mehg_ratio ~ mean_fork_length_mm, data = yp_df)
summary(fork_lm2)     # not significant

####################################################################
#### ROM meta-analysis and meta-regression with weight ####

yp_all_rom <- metacont(
  data      = yp_df,
  n.c       = thg_n,
  mean.c    = thg_mean_ug_g,
  sd.c      = thg_sd,
  n.e       = mehg_n,
  mean.e    = mehg_mean_ug_g,
  sd.e      = mehg_sd,
  sm        = "ROM",
  backtransf = TRUE)

yp_all_rom

yp.weight.reg <- metareg(yp_all_rom, ~ mean_weight_g )
yp.weight.reg
bubble(yp.weight.reg, studlab = FALSE, backtransf = TRUE)

yp.age.reg <- metareg(yp_all_rom, ~ mean_age )
yp.age.reg
bubble(yp.age.reg, studlab = FALSE, backtransf = TRUE)

yp.fork.reg <- metareg(yp_all_rom, ~ mean_fork_length_mm )
yp.fork.reg
bubble(yp.fork.reg, studlab = FALSE, backtransf = TRUE)

####################################################################
#### Combined figure: sex + covariates (Age, Weight, Fork length)####

composite_plot <- ggarrange(
  yp_boxplot, age_plot, weight_plot, fork_plot,
  labels      = c("A", "B", "C", "D"),
  vjust       = 1,
  align       = "v",
  ncol        = 2,
  nrow        = 2,
  font.label  = list(size = 20, color = "black", face = "bold"),
  common.legend = TRUE)

ggsave2(
  "composite_yellow_perch_covary.jpg",
  composite_plot,
  width = 15, height = 10, dpi = 300)
