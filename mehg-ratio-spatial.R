#####################################################################
# Assessing spatial patterns of MeHg:THg Ratios
# By: Kristin Eccles
# Updated: Dec 9th, 2025
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
library(forcats)

####################################################################
#load data
df_lat <- read.csv("thg_mehg_metanalysis.csv", header=TRUE, fileEncoding="latin1")
df_lat$factor <- paste(df_lat$common_name, df_lat$tissue, sep=" ")

df_lat$factor <- paste(df_lat$common_name, df_lat$tissue, sep=" ")
df_lat$latitude<- as.numeric(df_lat$latitude)
df_lat$longitude<- as.numeric(df_lat$longitude)

# Remove rows where either latitude OR longitude is NA
df_lat <- df_lat[!is.na(df_lat$latitude) & !is.na(df_lat$longitude), ]

df_lat <- df_lat%>%
  filter( ! tissue == "plant",
          ! tissue == "fat",
          ! tissue == "heart",
          ! tissue == "skin",
          ! tissue == "whole body")

tissue_colors <- c(
  "muscle"        = "#4477AA",  # muted blue
  "liver"         = "#44AA99",  # teal-green
  "brain"         = "#AA3377",  # muted magenta
  "kidney"        = "#88CCEE",  # light cyan
  "egg"           = "#EEAA33"  # amber
)

#####################################################################
#### Modification by location ####
All_rom <-metacont(data=df_lat, 
                           n.c = thg_n,
                           mean.c = thg_mean_ug_g,
                           sd.c = thg_sd,
                           n.e = mehg_n,
                           mean.e = mehg_mean_ug_g,
                           sd.e =  mehg_sd,
                           sm = "ROM",
                           backtransf = TRUE)
## Mean-center latitude
All_rom$data$latitude_c <- scale(All_rom$data$latitude, center = TRUE, scale = FALSE)


## Meta-regression with centered latitude and interaction
All.reg <- metareg(All_rom, ~ latitude_c * tissue)
All.reg

## Bubble plot
#bubble(All.reg, studlab = FALSE, backtransf = TRUE)

####################################################################
#### Polar Bear Muscle ####
df_lat_pb <- subset(df_lat, df_lat$factor == "Polar bear muscle")
df_lat_pb

mean(df_lat_pb$mehg_mean_ug_g / df_lat_pb$thg_mean_ug_g)

pb_all_rom <- metacont(
  data = df_lat_pb,
  n.c = thg_n,
  mean.c = thg_mean_ug_g,
  sd.c = thg_sd,
  n.e = mehg_n,
  mean.e = mehg_mean_ug_g,
  sd.e = mehg_sd,
  sm = "ROM",
  backtransf = TRUE
)
pb_all_rom

pb.reg <- metareg(pb_all_rom, ~latitude)
pb.reg

bubble(pb.reg, studlab = FALSE, backtransf = TRUE)

plot_pb <- data.frame(
  latitude = pb_all_rom$data$latitude,
  TE       = pb_all_rom$TE,      # log ROM
  weight   = pb_all_rom$w.random
) %>%
  mutate(
    ROM = exp(TE),
    weight_scaled = weight / max(weight, na.rm = TRUE) * 20
  )

# Regression line from metareg
pb_beta0 <- coef(pb.reg)[1]
pb_beta1 <- coef(pb.reg)[2]

reg_pb <- data.frame(
  latitude = seq(
    min(plot_pb$latitude, na.rm = TRUE),
    max(plot_pb$latitude, na.rm = TRUE),
    length.out = 100
  )
) %>%
  mutate(
    ROM = exp(pb_beta0 + pb_beta1 * latitude)
  )

pb_muscle_plot <- ggplot(plot_pb, aes(x = latitude, y = ROM)) +
  geom_point(
    aes(size = weight_scaled),
    shape = 21,
    fill = "steelblue",
    alpha = 0.6
  ) +
  geom_line(
    data = reg_pb,
    linewidth = 1,
    color = "black"
  ) +
  scale_size_identity() +
  labs(title = "Polar Bear Muscle",
    x = "Latitude (Degrees)",
    y = "Ratio of Means (MeHg/tHg)"
  ) +
  theme_classic(base_size = 14)

pb_muscle_plot

####################################################################
#### Walleye Muscle ####
df_lat_walleye <-subset(df_lat, df_lat$factor =="Walleye muscle")
df_lat_walleye
Walleye_all_rom <-metacont(data=df_lat_walleye, 
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

# plot
bubble(Walleye.reg, studlab = FALSE, backtransf = TRUE)

plot_walleye <- data.frame(
  latitude = Walleye_all_rom$data$latitude,
  TE       = Walleye_all_rom$TE,      # log ROM
  weight   = Walleye_all_rom$w.random
) %>%
  mutate(
    ROM = exp(TE),
    weight_scaled = weight / max(weight, na.rm = TRUE) * 20)

# Regression line from metareg
walleye_beta0 <- coef(Walleye.reg)[1]
walleye_beta1 <- coef(Walleye.reg)[2]

reg_walleye <- data.frame(
  latitude = seq(
    min(plot_walleye$latitude, na.rm = TRUE),
    max(plot_walleye$latitude, na.rm = TRUE),
    length.out = 100
  )
) %>%
  mutate(
    ROM = exp(walleye_beta0 + walleye_beta1 * latitude)
  )

walleye <- ggplot(plot_walleye, aes(x = latitude, y = ROM)) +
  geom_point(
    aes(size = weight_scaled),
    shape = 21,
    fill = "steelblue",
    alpha = 0.6) +
  geom_line(
    data = reg_walleye,
    linewidth = 1,
    color = "black") +
  scale_size_identity() +
  labs(title = "Walleye Muscle",
    x = "Latitude (Degrees)",
    y = "Ratio of Means (MeHg/tHg)") +
  theme_classic(base_size = 14)

walleye 
####################################################################
#### Beluga ####
df_lat_Beluga <-subset(df_lat, df_lat$factor =="Beluga whale muscle")
df_lat_Beluga
Beluga_all_rom <-metacont(data=df_lat_Beluga, 
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

# Publication Plot
plot_beluga <- data.frame(
  latitude = Beluga_all_rom$data$latitude,
  TE       = Beluga_all_rom$TE,        # log ROM
  weight   = Beluga_all_rom$w.random
) %>%
  mutate(
    ROM = exp(TE),
    weight_scaled = weight / max(weight, na.rm = TRUE) * 20
  )

# Regression line from metareg
beluga_beta0 <- coef(Beluga.reg)[1]
beluga_beta1 <- coef(Beluga.reg)[2]

reg_beluga <- data.frame(
  latitude = seq(
    min(plot_beluga$latitude, na.rm = TRUE),
    max(plot_beluga$latitude, na.rm = TRUE),
    length.out = 100)) %>%
  mutate(ROM = exp(beluga_beta0 + beluga_beta1 * latitude))

beluga_plot <- ggplot(plot_beluga, aes(x = latitude, y = ROM)) +
  geom_point(
    aes(size = weight_scaled),
    shape = 21,
    fill = "steelblue",
    alpha = 0.6) +
  geom_line(
    data = reg_beluga,
    linewidth = 1,
    color = "black") +
  scale_size_identity() +
  labs( title = "Beluga Whale Muscle",
    x = "Latitude (Degrees)",
    y = "Ratio of Means (MeHg/tHg)"
  ) +
  theme_classic(base_size = 14)

beluga_plot

####################################################################
#### Ringed Seal Liver ####
df_lat_Ring_Seal_liver <- subset(df_lat, df_lat$factor == "Ringed seal liver")
df_lat_Ring_Seal_liver

Ring_liver_all_rom <- metacont(
  data = df_lat_Ring_Seal_liver,
  n.c = thg_n,
  mean.c = thg_mean_ug_g,
  sd.c = thg_sd,
  n.e = mehg_n,
  mean.e = mehg_mean_ug_g,
  sd.e = mehg_sd,
  sm = "ROM",
  backtransf = TRUE)
Ring_liver_all_rom

Ring_liver.reg <- metareg(Ring_liver_all_rom, ~latitude)
Ring_liver.reg

bubble(Ring_liver.reg, studlab = FALSE)

# Publication Plot
plot_seal <- data.frame(
  latitude = Ring_liver_all_rom$data$latitude,
  TE       = Ring_liver_all_rom$TE,      # log ROM
  weight   = Ring_liver_all_rom$w.random) %>%
  mutate(
    ROM = exp(TE),
    weight_scaled = weight / max(weight, na.rm = TRUE) * 20)

# Regression line from metareg
ring_beta0 <- coef(Ring_liver.reg)[1]
ring_beta1 <- coef(Ring_liver.reg)[2]

reg_ring <- data.frame(
  latitude = seq(
    min(plot_seal$latitude, na.rm = TRUE),
    max(plot_seal$latitude, na.rm = TRUE),
    length.out = 100)) %>%
  mutate(
    ROM = exp(ring_beta0 + ring_beta1 * latitude))

ring_liver_plot <- ggplot(plot_seal, aes(x = latitude, y = ROM)) +
  geom_point(
    aes(size = weight_scaled),
    shape = 21,
    fill = "steelblue",
    alpha = 0.6) +
  geom_line(
    data = reg_ring,
    linewidth = 1,
    color = "black"
  ) +
  scale_size_identity() +
  labs(title = "Ringed Seal Liver",
       x = "Latitude (Degrees)",
    y = "Ratio of Means (MeHg/tHg)") +
  theme_classic(base_size = 14)

ring_liver_plot

##############################################################################

combo_meta_reg <- ggarrange(beluga_plot, ring_liver_plot, pb_muscle_plot,walleye,
  ncol          = 2,
  nrow          = 2, 
  labels        = c("B","C", "D", "E"),
  font.label    = list(size = 20, face = "bold"))

combo_meta_reg2 <- ggarrange(ratio_boxplot,combo_meta_reg,
                            ncol          = 1,
                            nrow          = 2, 
                            labels        = c("A"," "),
                            font.label    = list(size = 20, face = "bold"))

combo_meta_reg2

# Export as TIFF with 300 DPI
ggsave(
  filename = "Figure2.tiff",
  plot     = combo_meta_reg2,
  width    = 15,
  height   = 15,
  dpi      = 300)


