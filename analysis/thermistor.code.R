library(ggplot2)
library(rockchalk)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom)
library(dplyr)
library(ggpubr)

##Colour blind friendly palette ##
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

setwd("~/Documents/cottonwood_drydown2023/leaf.thermistors")
thermistor.data <- read.csv("all.leaf.thermistor.data.csv")

str(thermistor.data)
thermistor.data$tree.id <- factor(thermistor.data$tree.id)
thermistor.data$population <- factor(thermistor.data$population)
thermistor.data$leaf <- factor(thermistor.data$leaf)
thermistor.data$logger.id <- factor(thermistor.data$logger.id)
thermistor.data$date <- as.Date(thermistor.data$date, "%m/%d/%y")

thermistor.data$population <- factor(thermistor.data$population, levels = c("CCR", "NRV", "TSZ", "JLA"))

## subset afternoon temperature data between 2pm-6pm
pm.thermistor.data <- thermistor.data[which(thermistor.data$time == "14:00:00" | thermistor.data$time == "14:30:00" | 
                                              thermistor.data$time == "15:00:00" | thermistor.data$time == "15:30:00" | 
                                              thermistor.data$time == "16:00:00" | thermistor.data$time == "16:30:00" | 
                                              thermistor.data$time == "17:00:00" | thermistor.data$time == "17:30:00" | 
                                              thermistor.data$time == "18:00:00"),]

#######################################################################################################################
#######################################################################################################################
############ report means and summary stats for leaf temperature 

pm.thermistor.cleaned.data %>%
  filter(date > "2023-08-25" & date <"2023-08-29") %>%
  group_by(population) %>%
  summarize(mean_leafT = mean(temp, na.rm = TRUE),
            max_leafT = max(temp, na.rm = TRUE),
            min_leafT = min(temp, na.rm = TRUE))


############ report means and summary stats for Tleaf data [IN SPECIFIC DATE RANGE] 
pm.thermistor.cleaned.data %>%
  select(population,date,delta.t) %>%
  filter(date > "2023-08-25" & date <"2023-08-29") %>%
  #group_by(population) %>%
  summarize(mean_deltaT = mean(delta.t, na.rm = TRUE),
            max_deltaT = max(delta.t, na.rm = TRUE),
            min_deltaT = min(delta.t, na.rm = TRUE))

#######################################################################################################################
#######################################################################################################################


ggplot(pm.thermistor.cleaned.data, 
       aes(x = elevation, 
           y = temp, 
           colour = elevation, 
           shape = elevation)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "Source elevation (m)", y = "Tleaf (°C)") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  scale_colour_manual(values = cbbPalette)

tleaf.mixed.model <- lmer(temp ~ elevation + (1|leaf), data = pm.thermistor.cleaned.data)
anova(tleaf.mixed.model)
summary(tleaf.mixed.model)
plot(tleaf.mixed.model)
emmeans(tleaf.mixed.model, list(pairwise ~ elevation), adjust = "tukey")


ggplot(pm.thermistor.data, aes(x = date, y = temp, colour = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "Date", y = "Avg leaf temperature (°C)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  scale_colour_manual(values = cbbPalette) +
  ggtitle("Leaf temperature (2pm-6pm)") +
  geom_vline(xintercept = as.numeric(as.Date("0023-08-25")), linetype="dashed", colour = "red", size=1) +
  geom_vline(xintercept = as.numeric(as.Date("0023-08-29")), linetype="dashed", colour = "red", size=1)


#
######### Fig 2 panel A ##########
#
thermistor.linegraph <- ggplot(pm.thermistor.cleaned.data, aes(x = date, y = temp, colour = elevation, fill = elevation)) +
  stat_summary(geom = "line", 
               fun = mean, 
               linewidth = 0.75) +
  stat_summary(geom ='ribbon', 
               fun.data = mean_cl_normal, 
               fun.args = list(conf.int = 0.95), 
               alpha = 0.2, 
               colour = NA) +
  labs(x = "Date", 
       y = expression(T[leaf]*" (°C)"),
       colour = "Source elevation (m)", 
       fill = "Source elevation (m)") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1)) +
  border() +
  scale_colour_manual(values = c("red", "orange", "deepskyblue", "navy")) +
  scale_fill_manual(values = c("red", "orange", "deepskyblue", "navy")) +
  scale_x_date(date_breaks = "weeks", date_labels = "%d %b") +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-11")), linetype="dotted", colour = "gray", size = 0.65) +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-18")), linetype="dotdash", colour = "gray50", size = 0.65) +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-25")), linetype="dashed", colour = "black", size = 0.65) +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-29")), linetype="longdash", colour = "black", size = 0.65)
  
  
  geom_rect(aes(xmin = as.Date("2023-08-11"), xmax = as.Date("2023-08-18"), ymin = -Inf, ymax = Inf),
    fill = "lightblue", alpha = 0.01, linetype = "blank") +
  geom_rect(aes(xmin = as.Date("2023-08-18"), xmax = as.Date("2023-08-25"), ymin = -Inf, ymax = Inf),
    fill = "khaki", alpha = 0.01, linetype = "blank") +
  geom_rect(aes(xmin = as.Date("2023-08-25"), xmax = as.Date("2023-08-29"), ymin = -Inf, ymax = Inf),
    fill = "lightcoral", alpha = 0.01, linetype = "blank")
  
  #geom_vline(xintercept = as.numeric(as.Date("2023-08-11")), linetype="dashed", colour = "orange2", size = 0.5) +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-18")), linetype="dashed", colour = "orangered", size = 0.5) +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-25")), linetype="dashed", colour = "orangered4", size = 0.5) +
 geom_vline(xintercept = as.numeric(as.Date("2023-08-29")), linetype="dashed", colour = "orangered4", size = 0.5) + 
  layer(geom = "rect",
    mapping = aes(xmin = as.Date("2023-08-11"), xmax = as.Date("2023-08-18"), ymin = -Inf, ymax = Inf, fill = "orange2", alpha = 0.05),
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,
    params = list(color = NA)) +
  layer(geom = "rect",
    mapping = aes(xmin = as.Date("2023-08-18"), xmax = as.Date("2023-08-25"), ymin = -Inf, ymax = Inf, fill = "orangered", alpha = 0.05),
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,
    params = list(color = NA)) +
  layer(geom = "rect",
    mapping = aes(xmin = as.Date("2023-08-25"), xmax = as.Date("2023-08-29"), ymin = -Inf, ymax = Inf, fill = "orangered4", alpha = 0.05),
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,
    params = list(color = NA))
  


write.csv(pm.thermistor.data, 'pm.thermistor.data.csv', row.names = F)

setwd("~/Documents/cottonwood_drydown2023/leaf.thermistors")
pm.thermistor.cleaned.data <- read.csv("pm.leaf.thermistor.data.cleaned.csv")

pm.thermistor.cleaned.data$tree.id <- factor(pm.thermistor.cleaned.data$tree.id)
pm.thermistor.cleaned.data$population <- factor(pm.thermistor.cleaned.data$population)
pm.thermistor.cleaned.data$leaf <- factor(pm.thermistor.cleaned.data$leaf)
pm.thermistor.cleaned.data$logger.id <- factor(pm.thermistor.cleaned.data$logger.id)
pm.thermistor.cleaned.data$date <- as.Date(pm.thermistor.cleaned.data$date, "%m/%d/%y")

pm.thermistor.cleaned.data$population <- factor(pm.thermistor.cleaned.data$population, levels = c("CCR", "NRV", "TSZ", "JLA"))

# Add elevation column
pm.thermistor.cleaned.data <- pm.thermistor.cleaned.data %>%
  mutate(elevation = case_when(
    population == "CCR" ~ 72,
    population == "NRV" ~ 666,
    population == "TSZ" ~ 1219,
    population == "JLA" ~ 1521,
    TRUE ~ NA_real_))

pm.thermistor.cleaned.data$elevation <- factor(pm.thermistor.cleaned.data$elevation)

######### combined thermistor and weather station data for 2pm-6pm, 25th July - 25th September
thermistor.wstation.pm.data <- read.csv("pm.leaf.thermistor.data.airtemp.cleaned.csv")

str(thermistor.wstation.pm.data)
thermistor.wstation.pm.data$tree.id <- factor(thermistor.wstation.pm.data$tree.id)
thermistor.wstation.pm.data$population <- factor(thermistor.wstation.pm.data$population)
thermistor.wstation.pm.data$leaf <- factor(thermistor.wstation.pm.data$leaf)
thermistor.wstation.pm.data$logger.id <- factor(thermistor.wstation.pm.data$logger.id)
thermistor.wstation.pm.data$date <- as.Date(thermistor.wstation.pm.data$date, "%m/%d/%y")

thermistor.wstation.pm.data$population <- factor(thermistor.wstation.pm.data$population, levels = c("CCR", "NRV", "TSZ", "JLA"))

# Add elevation column
thermistor.wstation.pm.data <- thermistor.wstation.pm.data %>%
  mutate(elevation = case_when(
    population == "CCR" ~ 72,
    population == "NRV" ~ 666,
    population == "TSZ" ~ 1219,
    population == "JLA" ~ 1521,
    TRUE ~ NA_real_))

thermistor.wstation.pm.data$elevation <- factor(thermistor.wstation.pm.data$elevation)

#######################################################################################################################
############ report means and summary stats for delta T data 

thermistor.wstation.pm.data %>%
  group_by(population) %>%
  summarize(mean_deltaT = mean(delta.t, na.rm = TRUE),
            max_deltaT = max(delta.t, na.rm = TRUE),
            min_deltaT = min(delta.t, na.rm = TRUE))

############ report means and summary stats for delta T data [IN SPECIFIC DATE RANGE] ###################
thermistor.wstation.pm.data %>%
  select(population,date,delta.t) %>%
  filter(date > "2023-08-25" & date <"2023-08-29") %>%
  #group_by(population) %>%
  summarize(mean_deltaT = mean(delta.t, na.rm = TRUE),
            max_deltaT = max(delta.t, na.rm = TRUE),
            min_deltaT = min(delta.t, na.rm = TRUE))



#
######### Fig 2 panel B ##########
#
deltaT.linegraph <- ggplot(thermistor.wstation.pm.data, aes(x = date, y = delta.t, colour = elevation, fill = elevation)) +
  stat_summary(fun = mean, 
               geom = "line", 
               linewidth = 0.75) +
  stat_summary(geom ='ribbon', 
               fun.data = mean_cl_normal, 
               fun.args = list(conf.int = 0.95), 
               alpha = 0.2, 
               colour = NA) +
  labs(x = "Date", 
       y = expression(T[leaf]*" - "*T[air]*" (°C)"), 
       colour = "Source elevation (m)", 
       fill = "Source elevation (m)") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1), 
        legend.position = c(0.9, 0.15)) +
  border() +
  scale_colour_manual(values = c("red", "orange", "deepskyblue", "navy")) +
  scale_fill_manual(values = c("red", "orange", "deepskyblue", "navy")) +
  scale_x_date(date_breaks = "weeks", date_labels = "%d %b") +
  geom_hline(yintercept = 0, linetype="dashed", colour = "black", size=0.5) + 
  geom_vline(xintercept = as.numeric(as.Date("2023-08-11")), linetype="dotted", colour = "gray", size = 0.65) +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-18")), linetype="dotdash", colour = "gray50", size = 0.65) +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-25")), linetype="dashed", colour = "black", size = 0.65) +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-29")), linetype="longdash", colour = "black", size = 0.65)
  
  
  #geom_rect(aes(xmin = as.Date("2023-08-11"), xmax = as.Date("2023-08-18"), ymin = -Inf, ymax = Inf),
           # fill = "lightblue", alpha = 0.01, linetype = "blank") +
 # geom_rect(aes(xmin = as.Date("2023-08-18"), xmax = as.Date("2023-08-25"), ymin = -Inf, ymax = Inf),
           # fill = "khaki", alpha = 0.01, linetype = "blank") +
  #geom_rect(aes(xmin = as.Date("2023-08-25"), xmax = as.Date("2023-08-29"), ymin = -Inf, ymax = Inf),
           # fill = "lightcoral", alpha = 0.01, linetype = "blank")
  #geom_vline(xintercept = as.numeric(as.Date("2023-08-11")), linetype="dashed", colour = "orange2", size = 0.5) +
  #geom_vline(xintercept = as.numeric(as.Date("2023-08-18")), linetype="dashed", colour = "orangered", size = 0.5) +
  #geom_vline(xintercept = as.numeric(as.Date("2023-08-25")), linetype="dashed", colour = "red", size = 0.5) +
  #geom_vline(xintercept = as.numeric(as.Date("2023-08-29")), linetype="dashed", colour = "red", size = 0.5) +
  #geom_hline(yintercept = 0, linetype="dashed", colour = "grey", size=0.5)

  
ggarrange(thermistor.linegraph, deltaT.linegraph, ncol = 1, nrow = 2, common.legend = TRUE, labels = "AUTO")


airT.linegraph <- ggplot(thermistor.wstation.pm.data, aes(x = date, y = air.temp)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.75) +
  labs(x = "Date", y = "Air temperature (°C)") +
  theme_classic(base_size = 12) +
 #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  scale_colour_manual(values = cbbPalette) +
  scale_x_date(date_breaks = "weeks", date_labels = "%d %b") +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-25")), linetype="dashed", colour = "red", size=1) +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-29")), linetype="dashed", colour = "red", size=1)

rh.linegraph <- ggplot(thermistor.wstation.pm.data, aes(x = date, y = rh)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.75) +
  labs(x = "Date", y = "Relative humidity (%)") +
  theme_classic(base_size = 12) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  scale_colour_manual(values = cbbPalette) +
  scale_x_date(date_breaks = "weeks", date_labels = "%d %b") +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-25")), linetype="dashed", colour = "red", size=1) +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-29")), linetype="dashed", colour = "red", size=1)

vpd.linegraph <- ggplot(thermistor.wstation.pm.data, aes(x = date, y = air.vpd.kpa)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.75) +
  labs(x = "Date", y = "VPD (kpa)") +
  theme_classic(base_size = 12) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  scale_colour_manual(values = cbbPalette) +
  scale_x_date(date_breaks = "weeks", date_labels = "%d %b") +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-25")), linetype="dashed", colour = "red", size=1) +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-29")), linetype="dashed", colour = "red", size=1)

ggarrange(airT.linegraph, rh.linegraph, vpd.linegraph, ncol = 1, nrow = 3, common.legend = TRUE)


temp.vpd.obs <- thermistor.wstation.pm.data %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(mean.leafT = mean(leaf.temp, na.rm = TRUE), 
                   mean.airT = mean(air.temp, na.rm = TRUE), 
                   mean.rh = mean(rh, na.rm = TRUE),
                   mean.vpd = mean(vpd.kpa, na.rm = TRUE))

write.csv(temp.vpd.obs, 'temp.vpd.observations.csv', row.names = F)


############################################################################################################
############################################################################################################
########weather station####
############################################################################################################
############################################################################################################
setwd("~/Documents/cottonwood_drydown2023/weather.station")
weather.data <- read.csv("230926_DBG cottonwood weather_Data15.csv")

weather.data$date <- as.Date(weather.data$date, "%m/%d/%y")

## subset weather data by half hours
halfhour.weather.data <- weather.data[which(weather.data$minute == "0" | weather.data$minute == "30"),]


## subset afternoon temperature data between 2pm-6pm
pm.halfhour.weather.data <- halfhour.weather.data[which(halfhour.weather.data$hour == "14" | halfhour.weather.data$hour == "15" | 
                                                    halfhour.weather.data$hour == "16" | halfhour.weather.data$hour == "17" | 
                                                    halfhour.weather.data$hour == "18"),]


write.csv(pm.halfhour.weather.data, 'pm.halfhour.weather.data.csv', row.names = F)



yuma.weather <- read.csv("Yuma Max temps_1997-2023_01.csv")


#Dot plot of days over 45°C at Yuma airport from 1997-2023
yuma.plot <- ggplot(yuma.weather, aes(x = year, y = total.days)) + 
  geom_smooth(aes(group = 1), 
              method = 'lm', 
              se = T, 
              colour = "black", 
              size = 0.5) +
  geom_point(size = 4, 
             colour = "red3") +
  stat_poly_eq(use_label(c("R2", "P")), 
               size = 5) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        plot.title = element_text(hjust = 1, 
                                  size = 16)) +
  labs(x = "Year", 
       y = "Number of days ≥ 45°C",
       title = "Yuma Airport") +
  scale_x_continuous(breaks = seq(1997, 2023, by = 2)) +
  border()



# First, calculate the mean total days for each year
mean_days <- yuma.weather %>% 
  group_by(year) %>% 
  summarize(mean_total_days = mean(total.days))

# Plot the bar chart with trendline
ggplot(yuma.weather, aes(x = year, y = total.days)) + 
  geom_bar(stat = "identity", fill = "lightblue", color = "black") +
  geom_smooth(data = mean_days, aes(x = year, y = mean_total_days), 
              method = 'lm', se = FALSE, colour = "red", size = 1) +
  labs(x = "Year", y = "Number of days ≥ 45°C", title = "Total Days ≥ 45°C per Year") +
  theme_classic(base_size = 16) +
  border()


phx.airport.weather <- read.csv("Phoenix Temp means_01_bp.csv")
str(phx.airport.weather)
phx.airport.weather$date <- as.Date(phx.airport.weather$date, "%y/%m/%d")

summer.2023 <- phx.airport.weather_2023 %>% 
  group_by(day) %>% 
  summarize(mean.max.temp = mean(mean.max.temp))

ggplot(phx.airport.weather, aes(x = day, y = mean.max.temp)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5), size = 0) +
  stat_summary(fun.data = mean_cl_normal, 
               na.rm = TRUE, 
               geom = "errorbar", 
               position = position_dodge(width = 0.5),
               alpha = 0.25) +
  geom_line(stat = "summary", fun = mean, color = "black", size = 1, group = 1) +
  theme_classic(base_size = 16) +
  labs(x = "Day of year", y = "Maximum air temperature (°C)") +
  scale_x_continuous(breaks = seq(180, 280, by = 20)) +
  border()

# Filter the data to include only observations from the year 2023
phx.airport.weather_2023 <- phx.airport.weather %>%
  filter(year == 2023)

# Plot the data with the current geom_line and the additional geom_line for the year 2023
phx.airport.weather.plot <- ggplot(phx.airport.weather, aes(x = date, y = mean.max.temp)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5), size = 0) +
  stat_summary(fun.data = mean_cl_normal, 
               na.rm = TRUE, 
               geom = "errorbar", 
               position = position_dodge(width = 0.5),
               alpha = 0.25) +
  geom_line(stat = "summary", fun = mean, 
            aes(color = "30-year mean"), 
            size = 1, 
            group = 1) +  # Line connecting mean values for each day
  geom_line(data = phx.airport.weather_2023, 
            aes(group = 1, 
                color = "2023"), 
            size = 1) +  # Additional line connecting mean values for each day in 2023
  scale_color_manual(values = c("30-year mean" = "black", "2023" = "red"), 
                     labels = c("2023", "30-year mean")) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        plot.title = element_text(hjust = 1, 
                                  size = 16)) +
  scale_x_date(date_breaks = "weeks", 
               date_labels = "%d %b",
               limits = as.Date(c('2023-07-03','2023-09-25'))) +
  labs(x = "Date", 
       y = "Maximum air temperature (°C)", 
       title = "Phoenix International Airport") +
  border() +
  theme(legend.position = c(0.15, 0.2),
        legend.title = element_blank())

ggarrange(yuma.plot, phx.airport.weather.plot, nrow = 2, labels = "AUTO")



################################ relative humidity
rh.data <- read.csv("phx.airport.rh.csv")

rh.data$date <- as.Date(rh.data$date, "%m/%d/%y")

mean_rh <- rh.data %>% 
  group_by(date) %>% 
  summarize(mean_pm_rh = mean(rh))

write.csv(mean_rh, 'pm.rh.data.csv', row.names = F)



############################################################################################################
############################################################################################################
########night time leaf temperature####
############################################################################################################
############################################################################################################
setwd("~/Documents/cottonwood_drydown2023/leaf.thermistors")
thermistor.data <- read.csv("all.leaf.thermistor.data.csv")

str(thermistor.data)
thermistor.data$tree.id <- factor(thermistor.data$tree.id)
thermistor.data$population <- factor(thermistor.data$population)
thermistor.data$leaf <- factor(thermistor.data$leaf)
thermistor.data$logger.id <- factor(thermistor.data$logger.id)
thermistor.data$date <- as.Date(thermistor.data$date, "%m/%d/%y")

thermistor.data$population <- factor(thermistor.data$population, levels = c("CCR", "NRV", "TSZ", "JLA"))


## subset night-time leaf temp data pm-6pm
night.thermistor.data <- thermistor.data[which(thermistor.data$time == "20:00:00" | thermistor.data$time == "20:30:00" | 
                                              thermistor.data$time == "21:00:00" | thermistor.data$time == "21:30:00" | 
                                              thermistor.data$time == "22:00:00" | thermistor.data$time == "22:30:00" | 
                                              thermistor.data$time == "23:00:00" | thermistor.data$time == "23:30:00" | 
                                              thermistor.data$time == "00:00:00" | thermistor.data$time == "00:30:00" | 
                                              thermistor.data$time == "01:00:00" | thermistor.data$time == "01:30:00" | 
                                              thermistor.data$time == "02:00:00" | thermistor.data$time == "02:30:00" | 
                                              thermistor.data$time == "03:00:00" | thermistor.data$time == "03:30:00" | 
                                              thermistor.data$time == "04:00:00" | thermistor.data$time == "04:30:00" | 
                                              thermistor.data$time == "05:00:00"),]

write.csv(night.thermistor.data, 'nighttime_tleaf.csv', row.names = F)


ggplot(night.thermistor.data, aes(x = date, y = temp, colour = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "Date", y = "Avg leaf temperature (°C)") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  scale_colour_manual(values = cbbPalette) +
  ggtitle("Night-time Tleaf (8pm-5am)") +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-11")), linetype="dotted", colour = "orange", size = 0.65) +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-18")), linetype="dotdash", colour = "darkorange", size = 0.65) +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-25")), linetype="dashed", colour = "red3", size = 0.65) +
  geom_vline(xintercept = as.numeric(as.Date("2023-08-29")), linetype="longdash", colour = "red3", size = 0.65)


night.thermistor.data %>%
  group_by(population) %>%
  summarize(mean_nightleafT = mean(temp, na.rm = TRUE),
            max_nightleafT = max(temp, na.rm = TRUE),
            min_nightleafT = min(temp, na.rm = TRUE))







