library(ggplot2)
library(rockchalk)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom)
library(dplyr)

##Colour blind friendly palette ##
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

setwd("~/Documents/cottonwood_drydown2023/water.potential")
wpotential.data <- read.csv("dbg_cottonwood_waterpotential.csv")

str(wpotential.data)
wpotential.data$pot <- factor(wpotential.data$pot)
wpotential.data$interval <- factor(wpotential.data$interval)
wpotential.data$treatment <- factor(wpotential.data$treatment)
wpotential.data$population <- factor(wpotential.data$population)
wpotential.data$pop.geno <- factor(wpotential.data$pop.geno)
wpotential.data$chamber <- factor(wpotential.data$chamber)
wpotential.data$time.of.day <- factor(wpotential.data$time.of.day)
wpotential.data$meas.week <- factor(wpotential.data$meas.week)
wpotential.data$genotype <- factor(wpotential.data$genotype)

wpotential.data$population <- factor(wpotential.data$population, levels = c("CCR-COL", "NRV-NEW", "TSZ-SAN", "JLA-JAK"))
wpotential.data$treatment <- factor(wpotential.data$treatment, levels = c("predrought", "drought", "postdrought"))

# Add elevation column
wpotential.data <- wpotential.data %>%
  mutate(elevation = case_when(
    population == "CCR-COL" ~ 72,
    population == "NRV-NEW" ~ 666,
    population == "TSZ-SAN" ~ 1219,
    population == "JLA-JAK" ~ 1521,
    TRUE ~ NA_real_))

wpotential.data$elevation <- factor(wpotential.data$elevation)


## subset data by time of day
predawn.wp.data <- wpotential.data[which(wpotential.data$time.of.day == "predawn"),]
midday.wp.data <- wpotential.data[which(wpotential.data$time.of.day == "midday"),]


elevation_colours <- c("red", "orange", "deepskyblue", "navy", "white")

predawn_fill <- ifelse(wpotential.data$time.of.day == 'predawn', 
                       elevation_colours[match(wpotential.data$elevation, unique(wpotential.data$elevation))], 
                       "midday")

predawn_fill_colors <- ifelse(wpotential.data$time.of.day == 'predawn', 
                              elevation_colours[match(wpotential.data$elevation, unique(wpotential.data$elevation))], 
                              NA)

wp.yearweek.fig <- 
  ggplot(wpotential.data, aes(x = week.of.year,
                              y = water.potential,
                              group = interaction(week.of.year, elevation, time.of.day),
                              colour = elevation,
                              shape = elevation,
                              fill = ifelse(time.of.day == "predawn", elevation, 
                                            ifelse(is.na(time.of.day), NA, "transparent")))) +  
  stat_summary(fun = mean, position = position_dodge(width = 0.5), size = 1) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE,  
               geom = "errorbar",  
               size = 1,
               position = position_dodge(width = 0.5),
               alpha = 0.7) +
  labs(x = "Week of year",  
       y = "Ψ (MPa)",  
       colour = "Source elevation (m)") +
  theme_classic(base_size = 20) +
  ylim(-3, 0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  scale_colour_manual(values = elevation_colours) +
  scale_shape_manual(name = "Source elevation (m)", values = c(23,22,21,24,25)) +
  scale_x_continuous(breaks = seq(28, 38, by = 1)) +
  scale_fill_manual(values = elevation_colours, na.value = "transparent") +
  guides(fill = FALSE) +
  geom_vline(xintercept = 32.5, linetype="dotted", colour = "orange", size = 0.5) +
  geom_vline(xintercept = 33.5, linetype="dotdash", colour = "darkorange", size = 0.5) +
  geom_vline(xintercept = 34.5, linetype="dashed", colour = "red", size = 0.5) +
  geom_vline(xintercept = 35.5, linetype="longdash", colour = "red3", size = 0.5)


#original code that works for supp figure without filling points for time of day
ggplot(wpotential.data, aes(x = week.of.year,
                            y = water.potential,
                            group = interaction(week.of.year, elevation, time.of.day),
                            colour = elevation,
                            shape = elevation)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5), size = 1) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", 
               size = 1,
               position = position_dodge(width = 0.5),
               alpha = 0.7) +
  labs(x = "Week of year", 
       y = "Leaf water potential (MPa)", 
       colour = "Source elevation (m)") +
  theme_classic(base_size = 20) +
  ylim(-3, 0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  scale_colour_manual(values = elevation_colours) +
  scale_shape_manual(name = "Source elevation (m)", values = c(23,22,21,24)) +
  scale_x_continuous(breaks = seq(28, 38, by = 1)) +
  geom_vline(xintercept = 32.5, linetype="dotted", colour = "orange", size = 0.5) +
  geom_vline(xintercept = 33.5, linetype="dotdash", colour = "darkorange", size = 0.5) +
  geom_vline(xintercept = 34.5, linetype="dashed", colour = "red", size = 0.5) +
  geom_vline(xintercept = 35.5, linetype="longdash", colour = "red3", size = 0.5)
##




ggplot(wpotential.data, aes(x = week.of.year, 
                            y = water.potential, 
                            group = interaction(time.of.day, week.of.year, elevation), 
                            colour = elevation, 
                            shape = elevation,
                            fill = ifelse(time.of.day == "predawn", "predawn", NA)
                           )
       ) +
  stat_summary(fun = mean,
               aes(fill = ifelse(time.of.day == "predawn", "predawn", NA)),
               position = position_dodge(width = 0.5), 
               size = 1) +
  stat_summary(fun.data = mean_cl_normal, 
               na.rm = TRUE, 
               geom = "errorbar", 
               size = 1,
               position = position_dodge(width = 0.5),
               alpha = 0.7) +
  scale_colour_manual(values = elevation_colours) +
  scale_shape_manual(values = c(5, 0, 1, 2)) +
  scale_fill_manual(values = c("predawn" = "pink")) +
  scale_x_continuous(breaks = seq(28, 38, by = 1))

################################################################################################
#visualising
ggplot(midday.wp.data, aes(x = week.of.year, y = water.potential)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5), size = 0.7) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", 
               position = position_dodge(width = 0.5),
               alpha = 0.5) +
  labs(x = "Week of year", y = "Water potential (MPa)", 
       colour = "Source elevation (m)", 
       shape = "Source elevation (m)") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  scale_colour_manual(values = cbbPalette[1:4]) +
  scale_shape_manual(name = "Elevation (m)", values = c(23, 22, 21, 24)) +
  scale_x_continuous(breaks = seq(28, 38, by = 1)) +
  geom_vline(xintercept = 32.5, linetype="dotted", colour = "orange", size = 0.5) +
  geom_vline(xintercept = 33.5, linetype="dotdash", colour = "darkorange", size = 0.5) +
  geom_vline(xintercept = 34.5, linetype="dashed", colour = "red", size = 0.5) +
  geom_vline(xintercept = 35.5, linetype="longdash", colour = "red3", size = 0.5)


bxp <- ggboxplot(predawn.wp.data, x = "week.of.year", y = "water.potential", add = "point")


###### one-way repeated measure ANOVA ##### 

#outlier test
predawn.wp.data %>%
  group_by(week.of.year) %>%
  identify_outliers(water.potential)

#normality assumption
predawn.wp.data %>%
  group_by(week.of.year) %>%
  shapiro_test(water.potential)

ggqqplot(predawn.wp.data, "water.potential", facet.by = "week.of.year")

#repeat measures anova
res.aov <- anova_test(data = predawn.wp.data, dv = water.potential, wid = pot, within = week.of.year)
get_anova_table(res.aov)

#post-hoc test
pwc <- predawn.wp.data %>%
  pairwise_t_test(
    water.potential ~ week.of.year, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc

pwc <- pwc %>% add_xy_position(x = "week.of.year")

bxp + 
  stat_pvalue_manual(pwc) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc))

####### Two-way repeated measures anova

predawn.summary <- predawn.wp.data %>%
  group_by(population, week.of.year) %>%
  get_summary_stats(water.potential, type = "mean_sd")


bxp.2 <- ggboxplot(
  predawn.wp.data, x = "week.of.year", y = "water.potential",
  color = "population", palette = "jco")
bxp.2

pd.outliers <- predawn.wp.data %>%
  group_by(population, week.of.year) %>%
  identify_outliers(water.potential)

predawn.wp.data %>%
  group_by(population, week.of.year) %>%
  shapiro_test(water.potential)

ggqqplot(predawn.wp.data, "water.potential", ggtheme = theme_bw()) +
  facet_grid(week.of.year ~ population, labeller = "label_both")

res.aov.2 <- anova_test(
  data = predawn.wp.data, dv = water.potential, wid = pot,
  within = c(population, week.of.year))
get_anova_table(res.aov.2)

model <- aov(water.potential ~ population*week.of.year + Error(pot/(population*week.of.year)), data = predawn.wp.data)
summary(model)


### Stats

#Pooled water potential
wp.mmodel <- lmer(water.potential ~ population*week.of.year*time.of.day + (1|pot) + (1|chamber), data = wpotential.data)
anova(wp.mmodel)
summary(wp.mmodel)
plot(wp.mmodel)

#Predawn water potential
pd.wp.mmodel <- lmer(water.potential ~ population*week.of.year + (1|pot) + (1|chamber), data = predawn.wp.data)
anova(pd.wp.mmodel)
summary(pd.wp.mmodel)
plot(pd.wp.mmodel)

#Midday water potential
md.wp.mmodel <- lmer(water.potential ~ population*week.of.year + (1|pot) + (1|chamber), data = midday.wp.data)
anova(md.wp.mmodel)
summary(md.wp.mmodel)
plot(md.wp.mmodel)

md.wp.mmodel.week <- lmer(water.potential ~ week.of.year.cat + (1|pot) + (1|chamber), data = midday.wp.data)
anova(md.wp.mmodel.week)

emmeans(md.wp.mmodel.week, list(pairwise ~ week.of.year.cat), adjust = "tukey")

##### midday water potential across measuring weeks
ggplot(midday.wp.data, aes(x = treatment, y = water.potential, colour = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "Treatment", y = "Water potential (MPa)") +
  ggtitle("Midday leaf water potential") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  scale_colour_manual(values = cbbPalette)

##### predawn water potential across measuring weeks
ggplot(predawn.wp.data, aes(x = treatment, y = water.potential, colour = population)) + 
  stat_summary(fun = mean, position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", position = position_dodge(width = 0.5)) +
  labs(x = "Treatment", y = "Water potential (MPa)") +
  ggtitle("Predawn leaf water potential") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border() +
  scale_colour_manual(values = cbbPalette)


ggplot(wpotential.data, aes(x = meas.week, y = water.potential, colour = time.of.day)) +
  geom_boxplot() +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  border()

predawn.anova <- aov(water.potential ~ population*treatment, data = predawn.wp.data)
summary(predawn.anova)

midday.anova <- aov(water.potential ~ population*treatment, data = midday.wp.data)
summary(midday.anova)


midday.wp.data$week.of.year.cat <- factor(midday.wp.data$week.of.year, levels = c("29", "32", "33", "34", "35", "36", "38"))

midday.wp.summary <- midday.wp.data %>%
  group_by(week.of.year) %>%
  summarize(n = n(),
            mean_wp = mean(water.potential, na.rm = TRUE),
            sd_wp = sd(water.potential, na.rm = TRUE),
            se_wp = sd_wp/sqrt(n))







