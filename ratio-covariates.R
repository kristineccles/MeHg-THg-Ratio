#####################################################################
# Effect of co-variates on MeHg:THg ratios - Age, Sex, Fork Length
# By: Kristin Eccles
# Updated: August 18th, 2023
# Written in R Version 4.3.1
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

# load data
yp_df <- read.csv("batchelar_yellow_perch.csv")

#####################################################################
# EDA
cor1 <- cor(yp_df[,c("mehg_ratio", "mean_age", "mean_weight_g", "mean_fork_length_mm")])
corrplot(cor1, addCoef.col = "black", tl.col = "black")

# generate a meta analysis 
m.gen <- metagen(TE = mehg_ratio,
                 seTE = mehg_ratio_se,
                 studlab = first_author,
                 data = yp_df,
                 sm = "SMD",
                 fixed = FALSE,
                 random = TRUE,
                 method.tau = "REML")
#### Sex ####
# test the effect of Sex
yp_boxplot <- ggplot(data= yp_df, aes(y = mehg_ratio, x = sex))+
  geom_boxplot()+
  theme_bw()+
  labs(y = "MeHg Ratio (%)", x = "Sex")
yp_boxplot

# mean difference in a meta-analysis
sex_effect <- update.meta(m.gen, subgroup = sex, tau.common = TRUE)
sex_effect

# simialr results with simple T.test
t.test(data= yp_df, mehg_ratio~sex)
#no difference between male and females

#### Age ####
# Age
age_plot <- ggplot(data= yp_df, aes(y = mehg_ratio, x = mean_age, color = sex))+
  geom_point()+
  geom_smooth(method=lm, se=TRUE)+
  theme_bw()+
  labs(y = "MeHg Ratio (%)", x = "Mean Age (years)", color = "Sex")
age_plot

age_lm <- lm(data = yp_df,  mehg_ratio~mean_age+sex)
summary(age_lm)
Anova(age_lm, type = 3)

age_lm2 <- lm(data =yp_df,  mehg_ratio~mean_age)
summary(age_lm2)

#fit the ANCOVA model
ancova_age <- aov(mehg_ratio~mean_age+sex, data = yp_df)
summary(ancova_age)

# Weight
mean_weight_g_plot <- ggplot(data= yp_df, aes(y = mehg_ratio, x = (mean_weight_g), color = sex))+
  geom_point()+
  geom_smooth(method=lm, se=TRUE)+
  theme_bw()+
  labs(y = "MeHg Ratio (%)", x = "Mean Weight (g)", color = "Sex")
mean_weight_g_plot

weight_lm <- lm(data =yp_df,  mehg_ratio~mean_weight_g + sex)
summary(weight_lm)
Anova(weight_lm, type = 3)

weight_lm2 <- lm(data =yp_df,  mehg_ratio~mean_weight_g)
summary(weight_lm2) #not significant

yp_all_rom <-metacont(data=yp_df, 
                           n.c = thg_n,
                           mean.c = thg_mean_ug_g,
                           sd.c = thg_sd,
                           n.e = mehg_n,
                           mean.e = mehg_mean_ug_g,
                           sd.e =  mehg_sd,
                           sm = "ROM",
                           backtransf = TRUE)
yp_all_rom
yp.weight.reg <- metareg(yp_all_rom, ~mean_weight_g+sex)
yp.weight.reg
bubble(yp.weight.reg, studlab = FALSE, backtransf = TRUE)

#### Fork Length ####
# fork_length
mean_fork_length_plot <- ggplot(data= yp_df, aes(y = mehg_ratio, x = mean_fork_length_mm, color = sex))+
  geom_point()+
  geom_smooth(method=lm, se=TRUE)+
  labs(y = "MeHg Ratio (%)", x = "Mean Fork Length (mm)", color = "Sex")+
  theme_bw()
mean_fork_length_plot

fork_length_lm <- lm(data = yp_df,  mehg_ratio~mean_fork_length_mm+sex)
summary(fork_length_lm)

fork_length_lm2 <- lm(data =yp_df,  mehg_ratio~mean_fork_length_mm)
summary(fork_length_lm2) #not significant

#All
all_lm <- lm(data =yp_df,  mehg_ratio~mean_age+ sex  + mean_weight_g)
summary(all_lm)
Anova(all_lm, type = 3)

composite_plot=ggarrange(yp_boxplot, age_plot, mean_weight_g_plot, mean_fork_length_plot,
                         labels = c( "A", "B", "C", "D"),
                         vjust = 1,
                         align = "v",
                         #hjust = -0.5,
                         ncol = 2, nrow = 2,
                         font.label = list(size = 20, color = "black", face = "bold"),
                         common.legend = TRUE)
cowplot::ggsave2("composite_yellow_perch_covary.jpg",composite_plot, width = 40, height = 25, dpi = 300)

