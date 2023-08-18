#####################################################################
# Assessing spatial patterns of MeHg:THg Ratios
# By: Kristin Eccles
# Updated: August 18th, 2023
# Written in R Version 4.3.1
#####################################################################

# load libraries 
library(ggplot2)
library(ggpubr)
library(cowplot)
library(sjPlot)
library(corrplot)
library(stats)
library(lmtest)
library(car)
library(multcomp)

#load data
df <- read.csv("thg_mehg_metanalysis.csv", header=TRUE, fileEncoding="latin1")
df$factor <- paste(df$common_name, df$tissue, sep=" ")
#df$factor <- gsub(df$factor, pattern = " ", replacement ="")
#edit data
df$factor <- paste(df$common_name, df$tissue, sep=" ")
df$latitude<- as.numeric(df$latitude)
df$longitude<- as.numeric(df$longitude)

#####################################################################
#### Modification by location ####
unique(df$factor)

#### All Muscle ####
All_rom <-metacont(data=df, 
                           n.c = thg_n,
                           mean.c = thg_mean_ug_g,
                           sd.c = thg_sd,
                           n.e = mehg_n,
                           mean.e = mehg_mean_ug_g,
                           sd.e =  mehg_sd,
                           sm = "ROM",
                           backtransf = TRUE)
All_rom
All.reg <- metareg(All_rom, ~latitude)
All.reg
bubble(All.reg, studlab = FALSE, backtransf = TRUE)

#### Walleye Muscle ####
df_walleye <-subset(df, df$factor =="Walleye muscle")
df_walleye
Walleye_all_rom <-metacont(data=df_walleye, 
                   n.c = thg_n,
                   mean.c = thg_mean_ug_g,
                   sd.c = thg_sd,
                   n.e = mehg_n,
                   mean.e = mehg_mean_ug_g,
                   sd.e =  mehg_sd,
                   sm = "ROM",
                   backtransf = TRUE)
Walleye_all_rom
Walleye.reg <- metareg(Walleye_all_rom, ~latitude)
Walleye.reg
bubble(Walleye.reg, studlab = FALSE, backtransf = TRUE)

# without outlier
df_walleye_subset <- df_walleye[c(1:3,5:12),]
Walleye_all_rom2 <-metacont(data=df_walleye_subset, 
                           n.c = thg_n,
                           mean.c = thg_mean_ug_g,
                           sd.c = thg_sd,
                           n.e = mehg_n,
                           mean.e = mehg_mean_ug_g,
                           sd.e =  mehg_sd,
                           sm = "ROM",
                           backtransf = TRUE)
summary(Walleye_all_rom2)
Walleye.reg2 <- metareg(Walleye_all_rom2, ~latitude)
Walleye.reg2
bubble(Walleye.reg2, studlab = FALSE, backtransf = TRUE)

#### Beluga ####
df_Beluga <-subset(df, df$factor =="Beluga whale muscle")
df_Beluga
Beluga_all_rom <-metacont(data=df_Beluga, 
                           n.c = thg_n,
                           mean.c = thg_mean_ug_g,
                           sd.c = thg_sd,
                           n.e = mehg_n,
                           mean.e = mehg_mean_ug_g,
                           sd.e =  mehg_sd,
                           sm = "ROM",
                           backtransf = TRUE)
Beluga_all_rom
Beluga.reg <- metareg(Beluga_all_rom, ~latitude)
Beluga.reg
bubble(Beluga.reg, studlab = FALSE)

#### Ring Seal liver ####
df_Ring_Seal_liver <-subset(df, df$factor =="Ringed seal liver")
df_Ring_Seal_liver
Ring_liver_all_rom <-metacont(data=df_Ring_Seal_liver, 
                          n.c = thg_n,
                          mean.c = thg_mean_ug_g,
                          sd.c = thg_sd,
                          n.e = mehg_n,
                          mean.e = mehg_mean_ug_g,
                          sd.e =  mehg_sd,
                          sm = "ROM",
                          backtransf = TRUE)
Ring_liver_all_rom
Ring_liver.reg <- metareg(Ring_liver_all_rom, ~latitude)
Ring_liver.reg
bubble(Ring_liver.reg, studlab = FALSE)

