rm(list = ls())

library(tidyverse)
library(rstatix)

# get growth curve modelling functions
source("scripts_2023/tecan_analysis_functions/tecan_funcs.r")


#### load data ####
# data is mean gray value and variation in gray values through time

density_S_SE <- read_csv('data/S_SE_density_data.csv')
density_SM_SEM <- read_csv('data/SM_SEM_density_data.csv')


variation_S_SE <- read_csv('data/S_SE_variation_data.csv')
variation_SM_SEM <- read_csv('data/SM_SEM_variation_data.csv')

# associate treatment id with Dish number (from python script)

key_S_SE <- tribble(
  ~Plate, ~treatment_rep,
  'Dish 1', 'S_1',
  'Dish 2', 'S_2',
  'Dish 3', 'SE_2',
  'Dish 4', 'SE_3',
  'Dish 5', 'S_3',
  'Dish 6', 'SE_1'
)

key_SM_SEM <- tribble(
  ~Plate, ~treatment_rep,
  'Dish 1', 'SM_3',
  'Dish 2', 'SM_1',
  'Dish 3', 'SEM_2',
  'Dish 4', 'SEM_3',
  'Dish 5', 'SM_2',
  'Dish 6', 'SEM_1'
)


#### wrangling ####


# pivot the lawn density data into long form 
df_S_SE <- density_S_SE %>% 
  left_join(key_S_SE, by = 'Plate') %>% 
  relocate(treatment_rep) %>% 
  select(-Plate) %>% 
  pivot_longer(cols = !treatment_rep, values_to = 'mean_grey', names_to = 'hour')


df_SM_SEM <- density_SM_SEM %>% 
  left_join(key_SM_SEM, by = 'Plate') %>% 
  relocate(treatment_rep) %>% 
  select(-Plate) %>% 
  pivot_longer(cols = !treatment_rep, values_to = 'mean_grey', names_to = 'hour')

# combine the two data frames
df <- bind_rows(df_S_SE, df_SM_SEM) %>% 
  separate(treatment_rep, sep = '_', into = c('species', 'replicate')) %>% 
  mutate(hour = as.numeric(hour),
         replicate = as.numeric(replicate)) %>% 
  group_by(species, replicate) %>% 
  mutate(
    normalized_mean_grey = mean_grey - min(mean_grey)
  ) %>% 
  ungroup()

lawn_growthcurve <- df %>% 
  ggplot(aes(x = hour, y = normalized_mean_grey, group = factor(replicate))) + 
  geom_line()+
  theme_bw(12) +
  labs(y = 'Lawn density \n (mean grey value)', x = 'Time (hours)') +
  theme(legend.position = 'none') +
  facet_wrap(~factor(species, levels = c('S', 'SE', 'SM','SEM')))  +
  theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          legend.position = 'bottom',         # Bottom right corner
          legend.justification = c(0, 1),
          legend.text = element_text(color="black", size = 12),
          axis.text = element_text(color="black", size = 12), 
          strip.text = element_text(color="black", size = 12))

saveRDS(lawn_growthcurve, file = 'rds_plots/lawn_growthcurve.rdata')

# make a new column called well to run Chacon's growth rate code
# smooth the data so that the curves can be modeled more easily
df2 <- df %>% 
  unite(col= 'well', c(species, replicate), sep = '_',remove = F) %>%
  group_by(well) %>%
  arrange(hour) %>%
  mutate(normal_grey_smooth = rollmedian(normalized_mean_grey, 11, fill = NA))

# fit growth rate models to curves

log_lin_gr <- fit_all_loglinear(df2, measure = 'normal_grey_smooth', tries = 50)
logistic_gr <-  fit_all_logistic(df2, measure = 'normal_grey_smooth', tries = 50)
baranyi_gr <- fit_all_baranyi(df2, measure = 'normal_grey_smooth', tries = 50)
loglin_lm_gr <- fit_all_loglinear_lm(df2, measure = 'normal_grey_smooth', surround = 5)  

left_join(df2, baranyi_gr) %>%
  rowwise() %>%
  mutate(OD_pred = baranyi(hour,r, lag, ymax, y0)) %>%
  ungroup() %>%
  plot_plate(measure = "OD_pred")+ # here is plot plate
  geom_line(aes(y = normal_grey_smooth), size = 2, color = "black")+
  geom_line(size = 1, color = "blue") +
  scale_y_log10()

left_join(df2, logistic_gr) %>%
  rowwise() %>%
  mutate(OD_pred = logistic(hour,r, lag, K, y0)) %>%
  ungroup() %>%
  plot_plate(measure = "OD_pred")+ # here is plot plate
  geom_line(aes(y = normal_grey_smooth), size = 2, color = "black")+
  geom_line(size = 1, color = "blue") +
  scale_y_log10()


## seems like the baranyi function fits the data pretty well for growth rate

lawn_growth_rate <- baranyi_gr %>% 
  separate(well, into = c('Species', 'Replicate'), sep = '_', remove = FALSE) %>% 
  ggplot(aes(x = factor(Species,levels = c('S','SE','SM','SEM')), y = r))+
  geom_boxplot() +
  labs(x = 'Community', y = 'Growth rate (per hour)') +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = 'bottom',         # Bottom right corner
        legend.justification = c(0, 1),
        legend.text = element_text(color="black", size = 12),
        axis.text = element_text(color="black", size = 12), 
        strip.text = element_text(color="black", size = 12))

saveRDS(lawn_growth_rate, file = 'rds_plots/lawn_growth_rate.rdata')

# normality of groups -- all groups are normal
baranyi_gr %>% 
  separate(well, into = c('Species', 'Replicate'), sep = '_', remove = FALSE) %>% 
  group_by(Species) %>% 
  shapiro_test(r)

# check homogeneity of variance
# p val > 0.05
baranyi_gr %>% 
  separate(well, into = c('Species', 'Replicate'), sep = '_', remove = FALSE) %>% 
  levene_test(r ~ Species)


# anova test 

gr_anova<- baranyi_gr %>% 
  separate(well, into = c('Species', 'Replicate'), sep = '_', remove = FALSE) %>% 
  anova_test(r ~ Species)

# ANOVA Table (type II tests)
# 
# Effect DFn DFd      F        p p<.05   ges
# 1 Species   3   8 43.342 2.71e-05     * 0.942

gr_anova<- baranyi_gr %>% 
  separate(well, into = c('Species', 'Replicate'), sep = '_', remove = FALSE) %>% 
  anova_test(r ~ Species)

gr_tukey <- baranyi_gr %>% 
  separate(well, into = c('Species', 'Replicate'), sep = '_', remove = FALSE) %>% 
  tukey_hsd(r ~ Species)






