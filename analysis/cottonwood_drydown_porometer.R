######## 
#####analyse porometer data ##########

library(ggplot2)
library(rockchalk)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom)
library(dplyr)
library(ggpubr)
library(ggpmisc)

##Colour blind friendly palette ##
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

setwd("~/Documents/cottonwood_drydown2023/porometer")
porometer.data <- read.csv("porometer.all.days.csv")

str(porometer.data)
porometer.data$pot <- factor(porometer.data$pot)
porometer.data$leaf <- factor(porometer.data$leaf)
porometer.data$population <- factor(porometer.data$population)
porometer.data$genotype <- factor(porometer.data$genotype)
porometer.data$pop.geno <- factor(porometer.data$pop.geno)
porometer.data$rep <- factor(porometer.data$rep)
porometer.data$am.pm <- factor(porometer.data$am.pm)
porometer.data$meas.week <- factor(porometer.data$meas.week)
porometer.data$treatment <- factor(porometer.data$treatment)
porometer.data$treatment.levels <- factor(porometer.data$treatment.levels)
porometer.data$date <- as.Date(porometer.data$date, "%m/%d/%y")

porometer.data$population <- factor(porometer.data$population, levels = c("CCR-COL", "NRV-NEW", "TSZ-SAN", "JLA-JAK"))
porometer.data$treatment.levels <- factor(porometer.data$treatment.levels, levels = c("predrought", "drought.1", "drought.2", "drought.3", "postdrought"))
porometer.data$treatment <- factor(porometer.data$treatment, levels = c("predrought", "drought", "postdrought"))

#porometer.data <- porometer.data %>%
  #rename(tleaf = tair)

# add column for delta T
#porometer.data <- porometer.data %>%
 # mutate(delta.t = t.ref-t.leaf)


# subset am vs pm data
am.porometer.data <- porometer.data[which(porometer.data$am.pm=="am"),]
pm.porometer.data <- porometer.data[which(porometer.data$am.pm!="am"),]

pm.porometer.data$elevation <- factor(pm.porometer.data$elevation)

treatment.pm.porometer.data <- pm.porometer.data[which(pm.porometer.data$treatment.levels!="drought.1"),]
treatment.pm.porometer.data <- treatment.pm.porometer.data[which(treatment.pm.porometer.data$treatment.levels!="drought.2"),]

treatment.pm.porometer.data$treatment.levels <- factor(treatment.pm.porometer.data$treatment.levels, levels = c("predrought", "drought.3", "postdrought"))

treatment.pm.porometer.data$elevation <- factor(treatment.pm.porometer.data$elevation)

treatment.pm.porometer.data <- treatment.pm.porometer.data %>%
  mutate(treatments.1 = case_when(
    treatment.levels == "predrought" ~ "Pre-drought",
    treatment.levels == "drought.3" ~ "Drought",
    treatment.levels == "postdrought" ~ "Post-drought"))

treatment.pm.porometer.data$treatments.1 <- factor(treatment.pm.porometer.data$treatments.1, levels = c("Pre-drought", "Drought", "Post-drought"))

pm.porometer.data <- pm.porometer.data %>%
  mutate(treatment.rename = case_when(
    treatment == "predrought" ~ "Pre-drought",
    treatment == "drought" ~ "Drought",
    treatment == "postdrought" ~ "Post-drought"))

pm.porometer.data$treatment.rename <- factor(pm.porometer.data$treatment.rename, levels = c("Pre-drought", "Drought", "Post-drought"))


ggplot(porometer.data, aes(x = date, y = daily.max, group = 1)) + 
  geom_point(na.rm = TRUE) +
  geom_line() +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#############################################################################
##### gsw x treatment (afternoon)
gsw.pm.plot <- ggplot(pm.porometer.data, aes(x = treatment.rename, 
                                             y = gsw, group = elevation, 
                                             colour = elevation, 
                                             shape = elevation)) + 
  stat_summary(fun = mean, 
               position = position_dodge(width = 0.5), 
               size = 1.5) +
  stat_summary(fun.data = mean_cl_normal, 
               na.rm = TRUE, 
               geom = "errorbar", 
               size = 1,
               position = position_dodge(width = 0.5),
               alpha = 0.7) +
  labs(x = "Treatment", 
  y = expression(paste("Afternoon ", italic("g"[sw])~(mol~m^{-2}~s^{-1}))), 
  colour = "Source
elevation (m)") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1)) +
  border() +
  scale_colour_manual(values = c("red", "orange", "deepskyblue", "navy")) +
  scale_shape_manual(name = "Source
elevation (m)", values = c(23,22,21,24))

##########################################################################################################################################
###### conductance vs week of year ###### SUPPLEMENTARY FIGURE ######

gs.yearweek.fig <- ggplot(pm.porometer.data, aes(x = week.of.year, 
                                                 y = gsw, 
                                                 group = elevation, 
                                                 colour = elevation, 
                                                 shape = elevation)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5), size = 1) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", 
               size = 1,
               position = position_dodge(width = 0.5),
               alpha = 0.7) +
  labs(x = "Week of year", y = expression(italic("g"[sw])~(mol~m^{-2}~s^{-1})), colour = "Source elevation (m)") +
  #labs(x = "Date", y = expression("Stomatal conductance"~(mol~m^{-2}~s^{-1})), colour = "Elevation (m)") +
  theme_classic(base_size = 20) +
  ylim(0, 0.2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  scale_colour_manual(values = elevation_colours) +
  scale_shape_manual(name = "Source elevation (m)", values = c(23,22,21,24)) +
  scale_x_continuous(breaks = seq(28, 38, by = 1)) +
  geom_vline(xintercept = 32.5, linetype="dotted", colour = "orange", size = 0.5) +
  geom_vline(xintercept = 33.5, linetype="dotdash", colour = "darkorange", size = 0.5) +
  geom_vline(xintercept = 34.5, linetype="dashed", colour = "red", size = 0.5) +
  geom_vline(xintercept = 35.5, linetype="longdash", colour = "red3", size = 0.5)
  #geom_rect(aes(xmin = as.Date("2023-08-11"), xmax = as.Date("2023-08-18"), ymin = -Inf, ymax = Inf),
            #fill = "lightblue", alpha = 0.01, linetype = "blank", position = "identity") +
  #geom_rect(aes(xmin = as.Date("2023-08-18"), xmax = as.Date("2023-08-25"), ymin = -Inf, ymax = Inf),
            #fill = "khaki", alpha = 0.008, linetype = "blank", position = "identity") +
  #geom_rect(aes(xmin = as.Date("2023-08-25"), xmax = as.Date("2023-08-29"), ymin = -Inf, ymax = Inf),
            #fill = "lightcoral", alpha = 0.008, linetype = "blank", position = "identity")
##########################################################################################################################################
##### 2 panel conductance and water potential vs. week of year ##### SUPPLEMENTARY FIGURE FOR PAPER

ggarrange(gs.yearweek.fig, wp.yearweek.fig, ncol = 2, nrow = 1, common.legend = TRUE, labels = "AUTO")


### Stats

gsw.mmodel <- lmer(gsw ~ population*week.of.year + (1|pot), data = pm.porometer.data)
anova(gsw.mmodel)
summary(gsw.mmodel)
plot(gsw.mmodel)

gsw.mmodel.2 <- lmer(gsw ~ population*treatment + (1|pot), data = pm.porometer.data)
anova(gsw.mmodel.2)
summary(gsw.mmodel.2)
plot(gsw.mmodel.2)

##########################################################################################################################################

#########    dual y-axes test no. 1
ggplot(pm.porometer.data, aes(x = date)) + 
  stat_summary(aes(y = gsw, colour = elevation), fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(aes(y = gsw, colour = elevation), fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  stat_summary(aes(y = t.ref), fun = mean, 
               geom = "line", position = position_dodge(width = 0.5), linetype = "solid") +
  labs(x = "Date", y = expression("Stomatal conductance"~(mol~m^{-2}~s^{-1})), 
       colour = "Elevation (m)") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  scale_colour_manual(values = cbbPalette) +
  scale_y_continuous(name = expression("Stomatal conductance"~(mol~m^{-2}~s^{-1})),
                     sec.axis = sec_axis(~./100, name = "Air Temperature (°C)"))

#########    dual y-axes test no. 2
# Value used to transform the data
coeff <- 100

# A few constants
temperatureColor <- "#69b3a2"
priceColor <- rgb(0.2, 0.6, 0.9, 1)

ggplot(data, aes(x=day)) +
  geom_line( aes(y=temperature), size=2, color=temperatureColor) + 
  geom_line( aes(y=price / coeff), size=2, color=priceColor) +
  scale_y_continuous(name = "Temperature (Celsius °)", sec.axis = sec_axis(~.*coeff, name="Price ($)"))



ggplot(treatment.pm.porometer.data, aes(x = elevation, y = gsw)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "Elevation (m)", y = expression(italic("g"[s])~(mol~m^{-2}~s^{-1})), colour = "Elevation (m)") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  ggtitle("Stomatal conductance (3 - 6pm)")

##### delta.t x treatment (afternoon)
ggplot(treatment.pm.porometer.data, aes(x = treatment.levels, y = delta.t, colour = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "Treatment", y = "Tair - Tleaf") +
  ggtitle("Afternoon Delta T") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  scale_colour_manual(values = cbbPalette)


##### gsw x elevation (afternoon)
ggplot(treatment.pm.porometer.data, aes(x = elevation, y = gsw, colour = elevation)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "Elevation (m)", y = "Stomatal conductance") +
  ggtitle("Afternoon stomatal conductance") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() 

##### delta.t x treatment (afternoon)
ggplot(treatment.pm.porometer.data, aes(x = elevation, y = delta.t, colour = elevation)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "Elevation (m)", y = "Tair - Tleaf") +
  ggtitle("Afternoon Delta T") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border()


gsw.aov.pm.model <- aov(gsw ~ treatment.levels*population, data = treatment.pm.porometer.data)
summary(gsw.aov.pm.model)

############################################################
#### summary stats #####
pivot1 <- pm.porometer.data %>%
  group_by(treatment, population) %>%
  summarize(mean_gsw = mean(gsw, na.rm = TRUE),
            sd_gsw = sd(gsw, na.rm = TRUE),
            n = n(), 
            se_gsw = sd_gsw/sqrt(n),
            mean_deltaT = mean(delta.t, na.rm = TRUE),
            sd_deltaT = sd(delta.t, na.rm = TRUE),
            se_deltaT = sd_deltaT/sqrt(n))


pm.porometer.data %>%
  group_by(treatment, elevation) %>%
  summarize(mean_gsw = mean(gsw, na.rm = TRUE))

############################################################

##### elevation x tleaf
ggplot(porometer.data, aes(x = elevation, y = t.leaf)) + 
  stat_summary(fun = mean) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar") +
  facet_grid(~am.pm) +
  labs(x = "Elevation (m)", y = "Leaf temperature (C)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border()




##### treatment x delta.t (afternoon)
deltaT.treatment.pm <- ggplot(treatment.pm.porometer.data, aes(x = treatment.levels, y = delta.t, colour = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "treatment", y = "Tair - Tleaf (C)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  ggtitle("Delta T (pm)") #+
  facet_grid(~meas.week)

##### treatment x Gs (afternoon)
gs.treatment.pm <- ggplot(treatment.pm.porometer.data, aes(x = treatment.levels, y = gsw, colour = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "treatment", y = "Gsw") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  ggtitle("Stomatal conductance (pm)") #+
  facet_grid(~meas.week)

##### treatment x delta.t (morning)
deltaT.treatment.am <- ggplot(treatment.pm.porometer.data, aes(x = treatment.levels, y = delta.t, colour = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "treatment", y = "Tair - Tleaf (C)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  ggtitle("Delta T (am)") #+
  facet_grid(~meas.week)
  
##### treatment x Gs (morning)
gs.treatment.am <- ggplot(treatment.pm.porometer.data, aes(x = treatment.levels, y = gsw, colour = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "treatment", y = "Gsw") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  ggtitle("Stomatal conductance (am)") #+
  facet_grid(~meas.week)
    
ggarrange(gs.elevation.am, deltaT.elevation.am, gs.elevation.pm, deltaT.elevation.pm, ncol = 2, nrow = 2)




##### DPLYR pivot table summary
mean.max.tleaf.pm <- pm.porometer.data %>%
  group_by(treatment.levels, pot, pop.geno, population, elevation, time, date, meas.week) %>%
  summarise(mean_tleaf = mean(t.leaf), max_tleaf = max(t.leaf))




# Mean max tleaf (afternoon) across treatments
max.tleaf.pm.vs.treatment <- ggplot(mean.max.tleaf.pm, aes(x = treatment.levels, y = max_tleaf, colour = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "treatment", y = "Max Tleaf (°C)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  ylim(35, 50)

# Mean tleaf (afternoon) across treatments
mean.tleaf.pm.vs.treatment <- ggplot(mean.max.tleaf.pm, aes(x = treatment.levels, y = mean_tleaf, colour = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "treatment", y = "Mean Tleaf (°C)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  ylim(35, 50)

ggarrange(mean.tleaf.pm.vs.treatment, max.tleaf.pm.vs.treatment, ncol = 2, nrow = 1, common.legend = TRUE)



##### treatment x tleaf (afternoon)
ggplot(pm.porometer.data, aes(x = treatment, y = t.leaf, colour = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "treatment", y = "Tleaf (°C)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  ggtitle("Mean Leaf temp (pm)")



##### delta.t vs. vpd.kpa (afternoon)
ggplot(pm.porometer.data, aes(x = vpd.kpa, y = delta.t, colour = meas.week)) + 
  stat_summary(fun = mean) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar") +
  labs(x = "VPD (kpa)", y = "Tair - Tleaf (C)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  ggtitle("Delta T vs. VPD (pm)")


ggplot(pm.porometer.data, aes(x = meas.week, y = rh)) + 
  stat_summary(fun = mean) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar") +
  #labs(x = "Elevation (m)", y = "Delta T (C)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border()



pm.porometer.data$elevation <- factor(pm.porometer.data$elevation)


###### Air temp vs leaf temp graphs ######
##### afternoon air temp vs leaf temp (all populations)
ggplot(pm.porometer.data, aes(x = t.ref, y = t.leaf, colour = elevation, shape = elevation)) + 
  geom_point(size = 1.8, stroke = 1, alpha = 0.75) +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic(base_size = 14) +
  scale_shape_manual(values = c(0, 1, 2, 6)) +
  scale_colour_manual(values = cbbPalette) +
  #scale_color_manual(values = c("predrought" = "orange", "drought" = "red", "postdrought" = "blue")) +
  labs(x = "Air temperature (°C)", y = "Leaf temperature (°C)", 
       colour = "Elevation (m)", shape = "Elevation (m)") +
  facet_wrap(~treatment) +
  xlim(35,51) +
  ylim(30,51) +
  border()


##### afternoon air temp vs leaf temp CCR
ccr.tleaf.vs.tair <- ggplot(ccr.porometer.pm.data, aes(x = t.ref, y = t.leaf, colour = treatment)) + 
  geom_point(size = 1.5) +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic(base_size = 18) +
  labs(x = "Air temperature (°C)", y = "Leaf temperature (°C)", colour = "treatment", shape = "treatment") +
  xlim(35,52) +
  ylim(30,52) +
  ggtitle("Afternoon Tleaf vs. Tair CCR")

##### afternoon air temp vs leaf temp JLA
jla.tleaf.vs.tair <- ggplot(jla.porometer.pm.data, aes(x = t.ref, y = t.leaf, colour = meas.week)) + 
  geom_point(size = 1.5) +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic(base_size = 18) +
  labs(x = "Air temperature (°C)", y = "Leaf temperature (°C)", colour = "Week", shape = "Week") +
  xlim(35,52) +
  ylim(30,52) +
  ggtitle("Afternoon Tleaf vs. Tair JLA")

##### afternoon air temp vs leaf temp NRV
nrv.tleaf.vs.tair <- ggplot(nrv.porometer.pm.data, aes(x = t.ref, y = t.leaf, colour = meas.week)) + 
  geom_point(size = 1.5) +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic(base_size = 18) +
  labs(x = "Air temperature (°C)", y = "Leaf temperature (°C)", colour = "Week", shape = "Week") +
  xlim(35,52) +
  ylim(30,52) +
  ggtitle("Afternoon Tleaf vs. Tair NRV")

##### afternoon air temp vs leaf temp TSZ
tsz.tleaf.vs.tair <- ggplot(tsz.porometer.pm.data, aes(x = t.ref, y = t.leaf, colour = meas.week)) + 
  geom_point(size = 1.5) +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic(base_size = 18) +
  labs(x = "Air temperature (°C)", y = "Leaf temperature (°C)", colour = "Week", shape = "Week") +
  xlim(35,52) +
  ylim(30,52) +
  ggtitle("Afternoon Tleaf vs. Tair TSZ")

ggarrange(ccr.tleaf.vs.tair, nrv.tleaf.vs.tair, tsz.tleaf.vs.tair, jla.tleaf.vs.tair, ncol = 2, nrow = 2)




##### afternoon air temperature across weeks
pm.air.temp <- ggplot(pm.porometer.data, aes(x = meas.week, y = t.ref)) + 
  stat_summary(fun = mean) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar") +
  labs(x = "Week", y = "Afternoon air temp (°C)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border()

##### afternoon humidity across weeks
pm.humidity <- ggplot(pm.porometer.data, aes(x = meas.week, y = rh)) + 
  stat_summary(fun = mean) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar") +
  labs(x = "Week", y = "Afternoon rel. humidity (%)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border()

##### afternoon vpd across weeks
pm.vpd <- ggplot(pm.porometer.data, aes(x = meas.week, y = vpd.kpa)) + 
  stat_summary(fun = mean) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar") +
  labs(x = "Week", y = "Afternoon leaf VPD (kpa)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border()

ggarrange(pm.air.temp, pm.humidity, pm.vpd, ncol = 3)

## humidity x Gsw
ggplot(pm.porometer.data, aes(x = rh, y = gsw, colour = population, shape = population)) + 
  geom_point() +
  stat_poly_line(se = FALSE) +
  stat_poly_eq(use_label(c("eq", "R2")), label.x = 0.9) +
  #geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 18) +
  labs(x = "Relative Humidity (%)", y = "Stomatal conductance") +
  ggtitle("Relative humidity vs. Stomatal conductance")

## Leaf temperature x Gsw
ggplot(pm.porometer.data, aes(x = t.leaf, y = gsw)) + 
  geom_point() +
  #stat_poly_line(se = FALSE) +
  #stat_poly_eq(use_label(c("eq", "R2")), label.x = 0.9) +
  #geom_smooth(method = "lm", se = FALSE) +
  stat_smooth(method = lm, formula = y ~ poly(x, 2)) +
  theme_classic(base_size = 18) +
  labs(x = "Leaf temperature (C)", y = "Stomatal conductance") +
  ggtitle("Afternoon Leaf temperature vs. Stomatal conductance")


###### Population specific graphs #########
#####
## subset data by population
ccr.porometer.data <- porometer.data[which(porometer.data$population=="CCR-COL"),]
jla.porometer.data <- porometer.data[which(porometer.data$population=="JLA-JAK"),]
nrv.porometer.data <- porometer.data[which(porometer.data$population=="NRV-NEW"),]
tsz.porometer.data <- porometer.data[which(porometer.data$population=="TSZ-SAN"),]


ccr.porometer.pm.data <- pm.porometer.data[which(pm.porometer.data$population=="CCR-COL"),]
jla.porometer.pm.data <- pm.porometer.data[which(pm.porometer.data$population=="JLA-JAK"),]
nrv.porometer.pm.data <- pm.porometer.data[which(pm.porometer.data$population=="NRV-NEW"),]
tsz.porometer.pm.data <- pm.porometer.data[which(pm.porometer.data$population=="TSZ-SAN"),]


##### CCR-COL  #### 

###gsw CCR
ggplot(ccr.porometer.pm.data, aes(x = meas.week, y = delta.t, colour = pop.geno, shape = pop.geno)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  #stat_poly_line(se = FALSE) +
  #stat_poly_eq(use_label(c("eq", "R2")), label.x = 0.9) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 18) #+
  labs(x = "Week", y = "Stomatal conductance") +
  ggtitle("CCR-COL (pm): Stomatal conductance")

## gsw x genotype CCR
ggplot(ccr.porometer.pm.data, aes(x = reorder(geno.new, +gsw), y = gsw)) +
  stat_summary(fun = mean) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar") +
  labs(x = "Genotype", y = "Gsw") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  ggtitle("CCR-COL: Genotype stomatal conductance")

ccr.porometer.model <- lmer(gsw ~ genotype*am.pm + (1|rep) + (1|date), data = ccr.porometer.data)
anova(ccr.porometer.model)
summary(ccr.porometer.model)
plot(ccr.porometer.model)

##### JLA-JAK  ####

## vpd x gsw JLA
ggplot(jla.porometer.pm.data, aes(x = meas.week, y = t.leaf, colour = pop.geno, shape = pop.geno)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  #stat_poly_line(se = FALSE) +
  #stat_poly_eq(use_label(c("eq", "R2")), label.x = 0.9) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 18) #+
  labs(x = "VPD (kPa)", y = "Stomatal conductance") +
  ggtitle("JLA-JAK: VPD vs. Stomatal conductance")

## gsw x genotype JLA
ggplot(jla.porometer.data, aes(x = reorder(pop.geno, +gsw), y = gsw)) +
  stat_summary(fun = mean) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar") +
  facet_grid(~am.pm) +
  labs(x = "Genotype", y = "Gsw") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  ggtitle("JLA-JAK: Genotype stomatal conductance")


jla.porometer.model <- lmer(gsw ~ geno.new*am.pm + (1|rep) + (1|date), data = jla.porometer.data)
anova(jla.porometer.model)
summary(jla.porometer.model)
plot(jla.porometer.model)

##### NRV-NEW  ####

## vpd x gsw NRV
ggplot(nrv.porometer.pm.data, aes(x = meas.week, y = t.leaf, colour = pop.geno, shape = pop.geno)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  #stat_poly_line(se = FALSE) +
  #stat_poly_eq(use_label(c("eq", "R2")), label.x = 0.9) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 18) #+
  labs(x = "VPD (kPa)", y = "Stomatal conductance") +
  ggtitle("NRV-NEW: VPD vs. Stomatal conductance")

## gsw x genotype NRV
ggplot(nrv.porometer.data, aes(x = reorder(geno.new, +gsw), y = gsw)) +
  stat_summary(fun = mean) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar") +
  facet_grid(~am.pm) +
  labs(x = "Genotype", y = "Gsw") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() #+
  ggtitle("NRV-NEW: Genotype stomatal conductance")

nrv.porometer.model <- lmer(gsw ~ geno.new*am.pm + (1|rep) + (1|date), data = nrv.porometer.data)
anova(nrv.porometer.model)
summary(nrv.porometer.model)
plot(nrv.porometer.model)

##### TSZ-SAN  ####

## vpd x gsw TSZ
ggplot(tsz.porometer.pm.data, aes(x = meas.week, y = gsw, colour = pop.geno, shape = pop.geno)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  #stat_poly_line(se = FALSE) +
  #stat_poly_eq(use_label(c("eq", "R2")), label.x = 0.9) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 18) #+
  labs(x = "VPD (kPa)", y = "Stomatal conductance") +
  ggtitle("TSZ-SAN: VPD vs. Stomatal conductance")

## gsw x genotype TSZ
ggplot(tsz.porometer.data, aes(x = reorder(geno.new, +gsw), y = gsw)) +
  stat_summary(fun = mean) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar") +
  facet_grid(~am.pm) +
  labs(x = "Genotype", y = "Gsw") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  ggtitle("TSZ-SAN: Genotype stomatal conductance")

tsz.porometer.model <- lmer(gsw ~ geno.new*am.pm + (1|rep) + (1|date), data = tsz.porometer.data)
anova(tsz.porometer.model)
summary(tsz.porometer.model)
plot(tsz.porometer.model) 

##################################
##################################
##################################


## vpd x gsw
ggplot(pm.porometer.data, aes(x = vpd.kpa, y = gsw)) + 
  geom_point() +
  #stat_poly_line(se = FALSE) +
  #stat_poly_eq(use_label(c("eq", "R2")), label.x = 0.9) +
  #geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 18) +
  labs(x = "ArVPD (kPa)", y = "Stomatal conductance") +
  ggtitle("Afternoon VPD vs. Stomatal conductance")

## vpd x delta T (afternoon)
ggplot(pm.porometer.data, aes(x = vpd.kpa, y = delta.t, colour = meas.week, shape = meas.week)) + 
  geom_point() +
  stat_poly_line(se = FALSE) +
  stat_poly_eq(use_label(c("eq", "R2")), label.x = 0.9) +
  geom_smooth(method = "lm", se = FALSE, aes(group=1)) +
  theme_classic(base_size = 18) +
  labs(x = "VPD (kPa)", y = "Tair - Tleaf") +
  ggtitle("VPD vs. Delta T (afternoon)")

## vpd x gsw (afternoon)
ggplot(pm.porometer.data, aes(x = vpd.kpa, y = gsw)) + 
  geom_point() +
  #stat_poly_line(se = FALSE) +
  #stat_poly_eq(use_label(c("eq", "R2")), label.x = 0.9) +
  #geom_smooth(method = "lm", se = FALSE, aes(group=1)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1, se =FALSE, aes(group=1)) +
  theme_classic(base_size = 18) +
  labs(x = "VPD (kPa)", y = "Stomatal conductance") +
  ggtitle("VPD vs. Stomatal conductance (afternoon)") +
  ylim(0, 0.7)

## air temp x gsw
ggplot(pm.porometer.data, aes(x = t.ref, y = gsw, colour = treatment)) + 
  geom_point() +
  stat_poly_line(se = FALSE) +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  #geom_smooth(method = "lm", se = FALSE, aes(group=1)) +
  theme_classic(base_size = 18) +
  labs(x = "Arvo Air temp (°C)", y = "Stomatal conductance") +
  ggtitle("Arvo Air temp vs. Stomatal conductance") +
  xlim(35, 50)

## delta T x gsw afternoon
ggplot(pm.porometer.data, aes(x = delta.t, y = gsw, colour = population, shape = population)) + 
  geom_point() +
  stat_poly_line(se = FALSE) +
  stat_poly_eq(use_label(c("eq", "R2")), label.x = 0.9) +
  #geom_smooth(method = "lm") +
  theme_classic(base_size = 18) +
  labs(x = "Tair - Tleaf", y = "Stomatal conductance") +
  ggtitle("Delta T vs. Stomatal conductance (afternoon)")

## vpd x gsw morning
ggplot(am.porometer.data, aes(x = vpd.kpa, y = gsw, colour = population, shape = population)) + 
  geom_point() +
  stat_poly_line(se = FALSE) +
  stat_poly_eq(use_label(c("eq", "R2")), label.x = 0.9) +
  #geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 18) +
  labs(x = "VPD (kPa)", y = "Stomatal conductance") +
  ggtitle("VPD vs. Stomatal conductance (morning)")


# add column for gsw:vpd ratio
porometer.data <- porometer.data %>%
  mutate(gsw.vpd.ratio = gsw/vpd.kpa)

porometer.data <- porometer.data %>%
  mutate(gsw.vpd.ratio = gsw/vpd.kpa)

##### vpd:gsw ratio x elevation
ggplot(treatment.pm.porometer.data, aes(x = treatment.levels, y = gsw.vpd.ratio, colour = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "Treatment", y = "Gsw:VPD ratio") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggtitle("Gsw:VPD ratio (afternoon)") +
  border() +
  scale_colour_manual(values = cbbPalette)

anova_model <- aov(gsw.vpd.ratio ~ population*am.pm, data = porometer.data)
summary(anova_model)



gsw.ancova_model <- aov(gsw ~ am.pm + population*vpd.kpa, data = porometer.data)
Anova(gsw.ancova_model, type = "III") 

#define the post hoc comparisons to make
postHocs <- glht(gsw.ancova_model, linfct = mcp(population = "Tukey"))
summary(postHocs)


porometer.data.model <- lmer(gsw ~ population + (vpd.kpa | population), data = porometer.data)
anova(porometer.data.model)
summary(porometer.data.model)
coef(porometer.data.model)$population
plot(porometer.data.model) 

ggplot(porometer.data, aes(x = am.pm, y = gsw)) + geom_boxplot()


##################################################################
########################################################################################
###### TCRIT ####################################################################################
########################################################################################
##################################################################

setwd("~/Documents/cottonwood_drydown2023/tcrit")
tcrit.data <- read.csv("all_tcrit_drydown2023.csv")

str(tcrit.data)
tcrit.data$pot <- factor(tcrit.data$pot)
tcrit.data$leaf <- factor(tcrit.data$leaf)
tcrit.data$population <- factor(tcrit.data$population)
tcrit.data$pop.geno <- factor(tcrit.data$pop.geno)
tcrit.data$rep <- factor(tcrit.data$rep)
tcrit.data$meas.week <- factor(tcrit.data$meas.week)
tcrit.data$treatment <- factor(tcrit.data$treatment)
tcrit.data$treatment.levels <- factor(tcrit.data$treatment.levels)
tcrit.data$elevation <- factor(tcrit.data$elevation)

tcrit.data$population <- factor(tcrit.data$population, levels = c("CCR-COL", "NRV-NEW", "TSZ-SAN", "JLA-JAK"))
tcrit.data$treatment.levels <- factor(tcrit.data$treatment.levels, levels = c("predrought", "drought.1", "drought.2", "drought.3", "postdrought"))
tcrit.data$treatment <- factor(tcrit.data$treatment, levels = c("predrought", "drought", "postdrought"))


tcrit.data <- tcrit.data %>%
  mutate(treatments.1 = case_when(
    treatment == "predrought" ~ "Pre-drought",
    treatment == "drought" ~ "Drought",
    treatment == "postdrought" ~ "Post-drought"))

tcrit.data$treatments.1 <- factor(tcrit.data$treatments.1, levels = c("Pre-drought", "Drought", "Post-drought"))

#### Tcrit means plotted against treatment, populations
ggplot(tcrit.data, aes(x = treatment.levels, y = tcrit, colour = population, shape = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "Treatment", y = "Tcrit (°C)") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  scale_colour_manual(values = cbbPalette) #+
  geom_hline(yintercept = 50.88, linetype = "dashed", colour = "#56B4E9") +
  geom_hline(yintercept = 49.36, linetype = "dashed", colour = "#E69F00") +
  geom_hline(yintercept = 49.31, linetype = "dashed", colour = "black") +
  geom_hline(yintercept = 49.06, linetype = "dashed", colour = "#009E73") 
  
  
tcrit.treatment.plot <- ggplot(tcrit.data, aes(x = treatments.1, 
                                               y = tcrit, group = elevation, 
                                               colour = elevation, 
                                               shape = elevation)) +
    stat_summary(fun = mean, 
                 position = position_dodge(width = 0.5),
                 size = 1.5) +
    stat_summary(fun.data = mean_cl_normal, 
                 na.rm = TRUE, 
                 geom = "errorbar", 
                 size = 1,
                 position = position_dodge(width = 0.5),
                 alpha = 0.7) +
    labs(x = "Treatment", 
         y = expression("T"[crit]~"(°C)"), 
         colour = "Source elevation (m)") +
    theme_classic(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, 
                                     vjust = 1, 
                                     hjust = 1)) +
    border() +
  scale_colour_manual(values = c("red", "orange", "deepskyblue", "navy")) +
  scale_shape_manual(name = "Source elevation (m)", values = c(23,22,21,24))

#### Mean thermal safety margins plotted against treatment, populations
safety.margin.plot <- ggplot(tcrit.data, aes(x = treatments.1, 
                                             y = tsm, group = elevation, 
                                             colour = elevation, 
                                             shape = elevation)) + 
    stat_summary(fun = mean, 
                 position = position_dodge(width = 0.5),
                 size = 1.5) +
    stat_summary(fun.data = mean_cl_normal, 
                 na.rm = TRUE, 
                 geom = "errorbar", 
                 size = 1,
                 position = position_dodge(width = 0.5),
                 alpha = 0.7) +
    labs(x = "Treatment", 
         y = expression(T[crit]*" - "*maximum~T[leaf]*" (°C)"), 
         colour = "Source elevation (m)") +
    theme_classic(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, 
                                     vjust = 1, 
                                     hjust = 1)) +
    border() +
  scale_colour_manual(values = c("red", "orange", "deepskyblue", "navy")) +
    scale_shape_manual(name = "Source elevation (m)", values = c(23,22,21,24)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "red")
  


  ggplot(tcrit.data, aes(x = population, y = tcrit, colour = population)) + 
    stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
    stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
                 geom = "errorbar", position = position_dodge(width = 0.5)) +
    labs(x = "Population", y = "Tcrit (°C)") +
    theme_classic(base_size = 15) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    border() +
    facet_wrap(~treatment.levels) +
    scale_colour_manual(values = cbbPalette)  
  
predrought.tcrit.data <- tcrit.data[which(tcrit.data$treatment.levels == "predrought"),]
drought.tcrit.data <- tcrit.data[which(tcrit.data$treatment.levels == "drought.3"),]
postdrought.tcrit.data <- tcrit.data[which(tcrit.data$treatment.levels == "postdrought"),]

##### DPLYR pivot table summaries
mean.tcrit.data.genotypes <- tcrit.data %>%
  dplyr::group_by(treatment, population, pop.geno) %>%
  dplyr::summarise(mean.tcrit = mean(tcrit, na.rm = TRUE), mean.tsm = mean(tsm, na.rm = TRUE))

mean.tcrit.data.populations <- tcrit.data %>%
  group_by(treatment, population) %>%
  summarise(n = n(),
            mean.tcrit = mean(tcrit, na.rm = TRUE), 
            sd.tcrit = sd(tcrit, na.rm  = TRUE),
            se.tcrit = sd.tcrit/sqrt(n),
            mean.tsm = mean(tsm, na.rm = TRUE),
            sd.tsm = sd(tsm, na.rm  = TRUE),
            se.tsm = sd.tsm/sqrt(n))
            

mean.tcrit.data.treatment <- tcrit.data %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarise(mean.tcrit = mean(tcrit, na.rm = TRUE), 
                   mean.tsm = mean(tsm, na.rm = TRUE))

pivot1 <- pm.porometer.data %>%
  group_by(treatment, population) %>%
  summarize(mean_gsw = mean(gsw, na.rm = TRUE),
            sd_gsw = sd(gsw, na.rm = TRUE),
            n = n(), 
            se_gsw = sd_gsw/sqrt(n),
            mean_deltaT = mean(delta.t, na.rm = TRUE),
            sd_deltaT = sd(delta.t, na.rm = TRUE),
            se_deltaT = sd_deltaT/sqrt(n))

### Individual graphs for each treatment with TSM lines

#Predrought
predrought.tcrit.tsm.plot <- ggplot(predrought.tcrit.data, aes(x = population, y = tcrit, colour = population, shape = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "Population", y = "Tcrit (°C)") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  ggtitle("Pre-drought") +
  ylim(38, 51) +
  scale_colour_manual(values = cbbPalette) +
  geom_hline(yintercept = 46.99, linetype = "dashed", color = "#56B4E9", size = 1) +
  geom_hline(yintercept = 44.92, linetype = "dashed", color = "#E69F00", size = 1) +
  geom_hline(yintercept = 44.62, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = 46.33, linetype = "dashed", color = "#009E73", size = 1) 

#Drought
drought.tcrit.tsm.plot <- ggplot(drought.tcrit.data, aes(x = population, y = tcrit, colour = population, shape = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "Population", y = "Tcrit (°C)") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  ggtitle("Peak drought") +
  ylim(38, 51) +
  scale_colour_manual(values = cbbPalette) +
  geom_hline(yintercept = 50.88, linetype = "dashed", color = "#56B4E9", size = 1) + #TSZ-SAN
  geom_hline(yintercept = 49.36, linetype = "dashed", color = "#E69F00", size = 1) + #NRV-NEW
  geom_hline(yintercept = 49.31, linetype = "dashed", color = "black", size = 1) + #CCR-COL
  geom_hline(yintercept = 49.06, linetype = "dashed", color = "#009E73", size = 1) #JLA-JAK

#Post drought
postdrought.tcrit.tsm.plot <- ggplot(postdrought.tcrit.data, aes(x = population, y = tcrit, colour = population, shape = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "Population", y = "Tcrit (°C)") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  ggtitle("Post-drought") +
  ylim(38, 51) +
  scale_colour_manual(values = cbbPalette) +
  geom_hline(yintercept = 38.97, linetype = "dashed", color = "#56B4E9", size = 1) + #TSZ-SAN
  geom_hline(yintercept = 38.99, linetype = "dashed", color = "#E69F00", size = 1) + #NRV-NEW
  geom_hline(yintercept = 38.23, linetype = "dashed", color = "black", size = 1) + #CCR-COL
  geom_hline(yintercept = 40.03, linetype = "dashed", color = "#009E73", size = 1) #JLA-JAK

ggarrange(predrought.tcrit.tsm.plot, drought.tcrit.tsm.plot, postdrought.tcrit.tsm.plot, ncol = 3, common.legend = TRUE)


### tcrit stats ###  
tcrit.anova.model <- aov(tcrit ~ treatment, data = tcrit.data)
summary(tcrit.anova.model)
emmeans(tcrit.anova.model, list(pairwise ~ treatment), adjust = "tukey")


tcrit.anova.model <- aov(tcrit ~ population*treatment, data = tcrit.data)
summary(tcrit.anova.model)


tsm.anova.model <- aov(tsm ~ population*treatment, data = tcrit.data)
summary(tsm.anova.model)

tcrit.mixed.model <- lmer(tcrit ~ population*treatment + (1|pot), data = tcrit.data)
anova(tcrit.mixed.model)
summary(tcrit.mixed.model)
plot(tcrit.mixed.model)

tsm.mixed.model <- lmer(tsm ~ population*treatment + (1|pot), data = tcrit.data)
anova(tsm.mixed.model)
summary(tsm.mixed.model)
plot(tsm.mixed.model)



##################################################################
########################################################################################
###### TSM & HSM figure ####################################################################################
########################################################################################
##################################################################

setwd("~/Documents/cottonwood_drydown2023")
tsm.hsm.data <- read.csv("tsm.hsm.data.csv")

str(tsm.hsm.data)
tsm.hsm.data$population <- factor(tsm.hsm.data$population, levels = c("ccr", "nrv", "tsz", "jla"))
tsm.hsm.data$elevation <- factor(tsm.hsm.data$elevation, levels = c("72", "666", "1219", "1521"))


TSM.TLPHSM <- 
  ggplot(tsm.hsm.data, aes(x = tlp.hsm, 
                           y = tsm, 
                           colour = elevation, 
                           shape = elevation,
                           group = 1)) + 
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = tsm - tsm.se, 
                    ymax = tsm + tsm.se, 
                    colour = elevation), 
                alpha = 0.6) +
  geom_errorbar(aes(xmin = tlp.hsm - tlp.hsm.se, 
                    xmax = tlp.hsm + tlp.hsm.se, 
                    colour = elevation), 
                alpha = 0.6) +
  geom_smooth(aes(group = 1), 
              method = 'lm', 
              se = FALSE, 
              colour = "black", 
              size = 0.8) +
  stat_poly_eq(use_label(c("R2", "P")),
               size = 5) +
  theme_classic(base_size = 20) +
  xlim(-1.5, 0.5) +
  scale_colour_manual(values = c("red", "orange", "deepskyblue", "navy")) +
  scale_shape_manual(name = "Source elevation (m)", 
                     values = c(23,22,21,24)) +
  labs(x = expression(Ψ[md]*" - "*Ψ[TLP]*" (MPa)"), 
       y = expression(T[crit]*" - "*maximum~T[leaf]*" (°C)"),
       colour = "Source elevation (m)") +
       border() 


TSM.p88HSM <- 
  ggplot(tsm.hsm.data, aes(x = p88.hsm, 
                           y = tsm, 
                           colour = elevation, 
                           shape = elevation, 
                           group = 1)) + 
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = tsm - tsm.se, 
                    ymax = tsm + tsm.se, 
                    colour = elevation), 
                alpha = 0.6) +
  geom_errorbar(aes(xmin = p88.hsm - p88.hsm.se, 
                    xmax = p88.hsm + p88.hsm.se, 
                    colour = elevation), 
                alpha = 0.6) +
  geom_smooth(aes(group = 1), 
              method = 'lm', 
              se = FALSE, 
              colour = "black", 
              size = 0.8) +  # Use group = 1 for an overall trend line
  stat_poly_eq(use_label(c("R2", "P")),
               size = 5) +
  theme_classic(base_size = 20) +
  xlim(-1.5, 0.5) +
  scale_colour_manual(values = c("red", "orange", "deepskyblue", "navy")) +
  scale_shape_manual(name = "Source elevation (m)", 
                     values = c(23,22,21,24)) +
  labs(x = expression(Ψ[md]*" - "*Ψ[88]*" (MPa)"), 
       y = expression(T[crit]*" - "*maximum~T[leaf]*" (°C)"), 
       colour = "Source elevation (m)") +
  border()



ggarrange(TSM.TLPHSM, TSM.p88HSM, ncol = 1, nrow = 2, common.legend = TRUE, labels = "AUTO")



##################################################################
########################################################################################
###### Fv/Fm ####################################################################################
########################################################################################
##################################################################

setwd("~/Documents/cottonwood_drydown2023/FvFm")
FvFm.data <- read.csv("FvFm.drydown.alldata.csv")

str(FvFm.data)
FvFm.data$pot <- factor(FvFm.data$pot)
FvFm.data$leaf <- factor(FvFm.data$leaf)
FvFm.data$population <- factor(FvFm.data$population)
FvFm.data$genotype <- factor(FvFm.data$genotype)
FvFm.data$rep <- factor(FvFm.data$rep)
FvFm.data$treatment <- factor(FvFm.data$treatment)
FvFm.data$date <- as.Date(FvFm.data$date, "%m/%d/%y")

FvFm.data <- FvFm.data %>%
  mutate(elevation = case_when(
    population == "CCR" ~ 72,
    population == "NRV" ~ 666,
    population == "TSZ" ~ 1219,
    population == "JLA" ~ 1521,
    TRUE ~ NA_real_))

FvFm.data$elevation <- factor(FvFm.data$elevation)

FvFm.data$population <- factor(FvFm.data$population, levels = c("CCR", "NRV", "TSZ", "JLA"))
FvFm.data$treatment <- factor(FvFm.data$treatment, levels = c("predrought", "drought", "postdrought"))

FvFm.data <- FvFm.data %>%
  mutate(treatments.1 = case_when(
    treatment == "predrought" ~ "Pre-drought",
    treatment == "drought" ~ "Drought",
    treatment == "postdrought" ~ "Post-drought"))

FvFm.data$treatments.1 <- factor(FvFm.data$treatments.1, levels = c("Pre-drought", "Drought", "Post-drought"))


fvfm.plot <- ggplot(FvFm.data, aes(x = treatments.1, 
                                   y = fv.fm, group = elevation, 
                                   colour = elevation, 
                                   shape = elevation)) + 
  stat_summary(fun = mean, na.rm = TRUE, 
               position = position_dodge(width = 0.5), 
               size = 1.5) +
  stat_summary(fun.data = mean_cl_normal, 
               na.rm = TRUE, 
               geom = "errorbar", 
               position = position_dodge(width = 0.5),
               size = 1,
               alpha = 0.7) +
  labs(x = "Treatment", 
       y = expression("Photosystem II efficiency"~"(F"[v]~"/ F"[m]~")"),
       colour = "Source elevation (m)") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1), 
        legend.position = c(0.9, 0.15)) +
  border() +
  scale_colour_manual(values = c("red", "orange", "deepskyblue", "navy")) +
  scale_shape_manual(name = "Source elevation (m)", values = c(23,22,21,24))


### 3 panel graph - TSM, FvFM, Tcrit ###
ggarrange(tcrit.treatment.plot, safety.margin.plot, fvfm.plot, ncol = 3, common.legend = TRUE, labels = "AUTO")

### 4 panel graph - Gs, TSM, FvFM, Tcrit ###
ggarrange(gsw.pm.plot, tcrit.treatment.plot, safety.margin.plot, fvfm.plot, ncol = 4, common.legend = TRUE, legend = "top")

FvFm.tsm.data <- read.csv("FvFm.drydown.alldata.tsm.csv")

FvFm.tsm.data$pot <- factor(FvFm.tsm.data$pot)
FvFm.tsm.data$leaf <- factor(FvFm.tsm.data$leaf)
FvFm.tsm.data$population <- factor(FvFm.tsm.data$population)
FvFm.tsm.data$genotype <- factor(FvFm.tsm.data$genotype)
FvFm.tsm.data$rep <- factor(FvFm.tsm.data$rep)
FvFm.tsm.data$treatment <- factor(FvFm.tsm.data$treatment)
str(FvFm.tsm.data)

FvFm.tsm.drought <- FvFm.tsm.data[which(FvFm.tsm.data$treatment == "drought"),]


formula <- y ~ poly(x, 2, raw = TRUE)

ggplot(FvFm.tsm.data, aes(x = tsm, y = fv.fm, colour = treatment)) + 
  geom_point() +
  scale_colour_manual(values = c("red", "darkorange", "navy")) +
  stat_poly_line(formula = formula, 
                 se = TRUE, 
                 colour = "black") +
  stat_poly_eq(formula = formula,
               colour= "black",
               use_label(c("R2", "P")),
               label.y = "bottom", label.x = "right") +
  theme_classic(base_size = 14) +
  labs(x = "Thermal Safety Margin (°C)", 
       y = "Fv/Fm", 
       colour = "Treatment", 
       shape = "Treatment") +
  border()
#### No relationship between tsm and FvFm in any of the treatment windows ###


FvFm.data %>%
  group_by(treatment, elevation) %>%
  summarize(n = n(),
            mean_fv.fm = mean(fv.fm, na.rm = TRUE),
            sd.fv.fm = sd(fv.fm, na.rm  = TRUE),
            se.fv.fm = sd.fv.fm/sqrt(n),
            max_fv.fm = max(fv.fm, na.rm = TRUE),
            min_fv.fm = min(fv.fm, na.rm = TRUE))

mean.tcrit.data.populations <- tcrit.data %>%
  group_by(treatment, population) %>%
  summarise(n = n(),
            mean.tcrit = mean(tcrit, na.rm = TRUE), 
            sd.tcrit = sd(tcrit, na.rm  = TRUE),
            se.tcrit = sd.tcrit/sqrt(n),
            mean.tsm = mean(tsm, na.rm = TRUE),
            sd.tsm = sd(tsm, na.rm  = TRUE),
            se.tsm = sd.tsm/sqrt(n))


### fv/fm stats ###  
FvFm.model <- aov(fv.fm ~ population*treatment, data = FvFm.data)
summary(FvFm.model)

fv.fm.mixed.model <- lmer(fv.fm ~ population*treatment + (1|pot), data = FvFm.data)
anova(fv.fm.mixed.model)
summary(fv.fm.mixed.model)
plot(fv.fm.mixed.model)


›##### Combine Tcrit and porometer data

#subset porometer data down to pre
subset.pm.porometer.data <- pm.porometer.data[which(pm.porometer.data$meas.week == "1" | pm.porometer.data$meas.week == "5" | pm.porometer.data$meas.week == "7"),]

mean.tcrit.all.data <- tcrit.data %>%
  group_by(treatment.levels, pot, pop.geno, population, elevation, rep, date, meas.week) %>%
  summarise(mean_tcrit = mean(tcrit))

mean.porometer.data.subset <- subset.pm.porometer.data %>%
  group_by(treatment.levels, pot, pop.geno, population, elevation, rep, date, meas.week) %>%
  summarise(mean_tleaf = mean(t.leaf), max_tleaf = max(t.leaf))

write.csv(mean.porometer.data.subset, 'subset.porometer.data.csv', row.names = F)
write.csv(mean.tcrit.all.data, 'subset.tcrit.data.csv', row.names = F)

