rm(list = ls())

library(tidyverse)
library(rstatix)

# growth rate functions
source(
  "scripts_2023/tecan_analysis_functions/tecan_funcs.r"
)

# ending counts
counts <- read_csv('data/end_counts_ilvA_panC_20230417.csv')

#read map data from 96 well plate format
meta <-
  plater::read_plate(file = "data/plate_map_20230417.csv",
                     #convert 96 well plate map into tidy dataset
                     well_ids_column = "well",
                     sep = ",")

meta$well <- gsub("(?<![0-9])0+", "", meta$well, perl = TRUE) # take out padding 0 in well numbers (A01 --> A1) 


# read in the tecan data
OD <- read_tecan_csv("data/spent_assay_20230417_OD.csv", sep = ",", measure = "OD")

# subtracting the blanks

OD <-  OD %>% group_by(well) %>% 
  mutate(OD = OD - min(OD)) %>%
  mutate(OD = ifelse(OD <= 0.5/200, 0.5/200, OD)) %>%
  ungroup()

plot_plate(OD)+
  labs(title = "y starts around 0.5 / 200")


# join meta data with OD data 

OD_lab <- inner_join(OD, meta, "well")

OD_lab <- OD_lab %>% 
  separate(col = values, into = c("carbon", "spent_media_phase", "mutation" ), sep = "_", remove = F)

# smooth data with a median smoothing, window size 11
OD_lab = OD_lab %>%
  group_by(well) %>%
  arrange(hour) %>%
  mutate(OD_smooth = rollmedian(OD, 11, fill = NA))

# cap data at the maximum reached
OD_lab = OD_lab %>%
  group_by(well) %>%
  mutate(OD_capped = ifelse(cycle > cycle[OD == max(OD)][1], max(OD), OD))

plot_plate(OD_lab, measure = "OD_smooth")+labs(title = "smoothed OD")
plot_plate(OD_lab, measure = "OD_capped")+labs(title = "capped OD")


## growth rates with subtracted background blank that is smoothed and capped at max value
## 

loglin_gr = fit_all_loglinear(OD_lab, measure = "OD_smooth", tries = 50)

logistic_gr = fit_all_logistic(OD_lab, measure = "OD_smooth", tries = 50)

baranyi_gr <- fit_all_baranyi(OD_lab, measure = "OD_smooth", tries = 50)

loglin_lm_gr = fit_all_loglinear_lm(OD_lab, measure = "OD_smooth", surround = 5)


##
plot_plate_growthrates(loglin_gr, response = "r")
plot_plate_growthrates(baranyi_gr, response = "r")





left_join(OD_lab, logistic_gr) %>%
  rowwise() %>%
  mutate(OD_pred = logistic(hour, r, lag,K, y0)) %>%
  ungroup %>%
  plot_plate(measure = "OD_pred")+ # here is plot plate
  geom_line(aes(y = OD), size = 2, color = "black")+
  geom_line(size = 1, color = "blue")+
  scale_y_log10()


left_join(OD_lab, baranyi_gr) %>%
  rowwise() %>%
  mutate(OD_pred = baranyi(hour,r, lag, ymax, y0)) %>%
  ungroup %>%
  plot_plate(measure = "OD_pred")+ # here is plot plate
  geom_line(aes(y = OD), size = 2, color = "black")+
  geom_line(size = 1, color = "blue")+
  scale_y_log10()

# plot correlation between baranyi and loglin gr
loglin_gr %>%
  mutate(loglinear_gr = r) %>%
  select(well, loglinear_gr) %>%
  left_join(baranyi_gr %>%
              mutate(baranyi_gr = r) %>%
              select(well, baranyi_gr)) %>%
  ggplot(aes(x = loglinear_gr, y = baranyi_gr))+
  geom_point()+
  stat_smooth(method = "lm")+
  ylim(0,1)+
  xlim(0,0.4)



### Plotting 

# background removed blank data plotted:

OD_lab %>%
  filter(!is.na(mutation)) %>% 
  filter(!spent_media_phase %in% c('fresh', 'stat')) %>% 
  ggplot(aes(x = hour, y = OD , group = well, color = carbon))+
  geom_line()+
  labs(title = 'Midlog phase spent media')+
  facet_wrap(~ mutation) +
  scale_y_log10() +
  ylim(0.0001, 0.1) +
  theme_bw(12)






OD_lab %>%
  filter(!is.na(mutation)) %>% 
  filter(!spent_media_phase %in% c('fresh', 'mid')) %>% 
  ggplot(aes(x = hour, y = OD , group = well, color = carbon))+
  geom_line()+
  labs(title = 'Stationary phase spent media')+
  facet_grid(~ mutation )

OD_lab %>%
  filter(!is.na(mutation)) %>% 
  filter(!spent_media_phase %in% c('fresh')) %>% 
  ggplot(aes(x = hour, y = OD , group = well, color = carbon))+
  geom_line()+
  facet_grid(spent_media_phase~ mutation ) +
  scale_y_log10()

# plot yields 
OD_yield <- get_yields(OD_lab) %>% 
  left_join(meta, by = 'well') %>% 
  separate(col = values, into = c("carbon", "spent_media_phase", "mutation" ), sep = "_", remove = F) %>% 
  group_by(carbon, mutation, spent_media_phase) %>% 
  summarise(
    mean_yield = mean(yield),
    sd_yield = sd(yield)
  ) %>% 
  mutate(
    carbon = case_when(
      carbon == 'gal' ~ 'Filtered S galactose spent media',
      carbon == 'succ' ~ 'Filtered S succinate spent media'),
    mutation = case_when(
      mutation == 'ilvA' ~ 'S∆*ilvA*',
      mutation == 'panC' ~ 'S∆*panC*'
    )
    )

OD_yield %>% 
  filter(!is.na(mutation)) %>% 
  filter(!spent_media_phase %in% c('fresh', 'stat')) %>% 
  ggplot(aes(x = carbon, y = mean_yield, fill = carbon)) +
  geom_bar(stat = 'identity', color = 'black', position = position_dodge())+
  geom_errorbar(aes(ymin = mean_yield, ymax = mean_yield + sd_yield), width = 0.5)+
  labs(
    x = '', y = 'Max OD600'
  )+
  facet_wrap(~mutation) +
  theme_bw(12) +
  theme(legend.position = 'none') 

S_spent_media_ilvA_panC<- OD_yield %>% 
  filter(!is.na(mutation)) %>% 
  filter(!spent_media_phase %in% c('fresh', 'stat')) %>% 
  ggplot(aes(x = mutation, y = mean_yield, fill = mutation)) +
  geom_bar(stat = 'identity', color = 'black', position = position_dodge())+
  geom_errorbar(aes(ymin = mean_yield, ymax = mean_yield + sd_yield), width = 0.5)+
  labs(
    x = '', y = "Absorbance (OD~600~)"
  )+
  facet_wrap(~factor(carbon,
                     levels = c("Filtered S galactose spent media","Filtered S succinate spent media" ))) +
  theme_bw(12) +
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"), 
        axis.text.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown(),
        axis.text = element_text(color="black")) +
  scale_fill_manual(values = c('S∆*ilvA*' = '#00BA38', 'S∆*panC*' = '#619cff'))

saveRDS(object =S_spent_media_ilvA_panC, file = 'rds_plots/S_spent_media_ilvA_panC.rdata')

ggsave('S_spent_media_ilvA_panC.pdf',width = 5, height = 3)
####

cfu_yield <- counts %>% 
  left_join(meta, by = 'well') %>% 
  separate(col = values, into = c("carbon", "spent_media_phase", "mutation" ), sep = "_", remove = F) %>% 
  mutate(total_cfu = (count / (dilution * 0.005)) * .2) %>% # total cfu calculated by getting cfu/mL then multiplying the well volume (i.e., 0.2mL)
  group_by(carbon, mutation, spent_media_phase) %>% 
  summarise(
    mean_yield = mean(total_cfu),
    sd_yield = sd(total_cfu)
  ) %>% 
  mutate(
    carbon = case_when(
      carbon == 'gal' ~ 'Filtered S galactose spent media',
      carbon == 'succ' ~ 'Filtered S succinate spent media'),
    mutation = case_when(
      mutation == 'ilvA' ~ 'S∆*ilvA*',
      mutation == 'panC' ~ 'S∆*panC*'
    )
  )



cfu_yield %>% 
  filter(!is.na(mutation)) %>% 
  filter(!spent_media_phase %in% c('fresh', 'stat')) %>% 
  ggplot(aes(x = mutation, y = mean_yield, fill = mutation)) +
  geom_bar(stat = 'identity', color = 'black', position = position_dodge())+
  geom_errorbar(aes(ymin = mean_yield, ymax = mean_yield + sd_yield), width = 0.5)+
  labs(
    x = '', y = 'Total CFU'
  )+
  facet_wrap(~factor(carbon,
                     levels = c("Filtered S galactose spent media","Filtered S succinate spent media" ))) +
  theme_bw(12) +
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"), 
        axis.text.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown(),
        axis.text = element_text(color="black")) +
    scale_fill_manual(values = c('S∆*ilvA*' = '#00BA38', 'S∆*panC*' = '#619cff')) +
  scale_y_log10() 

