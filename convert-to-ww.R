#####################################################################
# Converting Database to wet weight (ww)
# By: Kristin Eccles
# Updated: August 18th, 2023
# Written in R Version 4.3.1
#####################################################################

# Load Libraries
library(ggplot2)
library(boot)
library(lmtest)
library(car)
library(arsenal)
library(dplyr)
library(tidyr)

# Load Data
df <- read.csv("tf_hg_mehg_database_clean.csv", header=TRUE, fileEncoding="latin1")
df_converted <- df

#####################################################################
# Prepare the data for analysis

# Convert SE to SD
# SE=SD/sqrt(n)
# SD = SE*sqrt(n)
df_converted <- mutate(df, c_thg_sd = ifelse(is.na(thg_sd), 
                                             (thg_se*sqrt(thg_n)), 
                                             thg_sd))
df_converted <- mutate(df_converted, c_thg_se = ifelse(is.na(thg_se), 
                                                       (thg_sd/sqrt(thg_n)), 
                                                       thg_se))

df_converted <- mutate(df_converted, c_mehg_sd = ifelse(is.na(mehg_sd), 
                                                        mehg_se*sqrt(mehg_n), 
                                                        mehg_sd))
df_converted <- mutate(df_converted, c_mehg_se = ifelse(is.na(mehg_se), 
                                                        mehg_sd/sqrt(mehg_n), 
                                                        mehg_se))

# Convert all the dry weight results to wet weight  
# Summarize moisture content by tissue and species
df_converted$group <- interaction(factor(df_converted$tissue),factor(df_converted$common_name), drop = T)
unique(df_converted$group)

tissue_species_summary <- aggregate(df_converted$moist_cont_perc, by=list(df_converted$group) , FUN= mean, na.rm=TRUE, na.action="na.pass")
colnames(tissue_species_summary) <-cbind("group", "avg_group_moisture")

#just by tissue
tissue_summary <- aggregate(df_converted$moist_cont_perc, by=list(df_converted$tissue) , FUN= mean, na.rm=TRUE, na.action="na.pass")
colnames(tissue_summary) <-cbind("tissue", "avg_tissue_moisture")

# add the average percentage to df
df_converted <- left_join(df_converted, tissue_species_summary, keep=FALSE)
df_converted <- left_join(df_converted, tissue_summary, keep=FALSE)

df_converted<- mutate(df_converted, c_moisture_percent  = ifelse(is.na(moist_cont_perc), 
                                                                 avg_group_moisture, 
                                                                 moist_cont_perc))
df_converted<- mutate(df_converted, c_moisture_percent  = ifelse(is.na(c_moisture_percent), 
                                                                 avg_tissue_moisture, 
                                                                 c_moisture_percent))


# Convert the concentrations with moisture content from dw to ww
# WW = ((100 - % of water)/100)* DW
df_converted <- mutate(df_converted, ww_thg_mean_ug_g = ifelse(wet_dry_wt=="dry", 
                                                               ((100 - c_moisture_percent)/100)* thg_mean_ug_g, 
                                                               thg_mean_ug_g))

df_converted <- mutate(df_converted, ww_thg_se = ifelse(wet_dry_wt=="dry", 
                                                        ((100 - c_moisture_percent)/100)* c_thg_se, 
                                                        c_thg_se))
df_converted <- mutate(df_converted, ww_thg_sd = ifelse(wet_dry_wt=="dry", 
                                                        ((100 - c_moisture_percent)/100)* c_thg_sd, 
                                                        c_thg_sd))

df_converted <- mutate(df_converted, ww_mehg_mean_ug_g = ifelse(wet_dry_wt=="dry", 
                                                                ((100 - c_moisture_percent)/100)* mehg_mean_ug_g, 
                                                                mehg_mean_ug_g))

df_converted <- mutate(df_converted, ww_mehg_se = ifelse(wet_dry_wt=="dry", 
                                                         ((100 - c_moisture_percent)/100)* c_mehg_se, 
                                                         c_mehg_se))
df_converted <- mutate(df_converted, ww_mehg_sd = ifelse(wet_dry_wt=="dry", 
                                                         ((100 - c_moisture_percent)/100)* c_mehg_sd, 
                                                         c_mehg_sd))

#Calculate ratio variables for analysis
write.csv(df_converted, "df_converted.csv")

#### Summary Table ####
summary(df_converted)
converted_hg_summary <- as.data.frame(df_converted) %>%
  group_by(df$common_name, df$tissue) %>%
  summarise(mean_thg = mean(ww_thg_mean_ug_g), 
            sd_thg=mean(ww_thg_sd, na.rm=T),
            sum_thg=sum(thg_n, na.rm=T), 
            mean_mehg = mean(ww_mehg_mean_ug_g, na.rm=T), 
            sd_mehg=mean(ww_mehg_sd, na.rm=T), 
            sum_mehg=sum(mehg_n, na.rm=T), 
            n=n())
converted_hg_summary <-as.data.frame(converted_hg_summary)
write.csv(converted_hg_summary, "converted_hg_summary.csv")
