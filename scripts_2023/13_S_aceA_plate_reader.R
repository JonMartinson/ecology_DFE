rm(list = ls())

library(tidyverse)

source("scripts_2023/tecan_analysis_functions/tecan_funcs.r")

#read map data from 96 well plate format
meta <- plater::read_plate(file = "data/aceA_retest_tecan_20230204_map.csv", #convert 96 well plate map into tidy dataset
                           well_ids_column = "well",
                           sep = ",")

meta$well <- gsub("(?<![0-9])0+", "", meta$well, perl = TRUE) # take out padding 0 in well numbers (A01 --> A1) 

# read in the tecan data

OD <- read_tecan_csv("data/aceA_full_OD_2030204.csv", sep = ",", measure = "OD")


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
    separate(col = values, into = c("mutation", "media_type"), sep = "_", remove = F)

# smooth data with a median smoothing, window size 11
OD_lab = OD_lab %>%
    group_by(well) %>%
    arrange(hour) %>%
    mutate(OD_smooth = rollmedian(OD, 11, fill = NA))

# cap data at the maximum reached
OD_lab = OD_lab %>%
    group_by(well) %>%
    mutate(OD_capped = ifelse(cycle > cycle[OD == max(OD)][1], max(OD), OD))

## growth rates with subtracted background blank that is smoothed and capped at max value
## 

#log linear growth rate model curve fitting
loglin_gr = fit_all_loglinear(OD_lab, measure = "OD_smooth", tries = 50)

### Plotting 

# background removed blank data plotted:

OD_lab2 <- OD_lab %>% 
    mutate(
        mutation = case_when(mutation == 'WT' ~ 'S WT',
                             mutation == 'aceA' ~ 'S∆*aceA*',
                             mutation == 'blank' ~ 'blank'),
        media_type = case_when(media_type == 'ace' ~ 'Acetate',
                          media_type == 'cocult' ~ 'Coculture with E',
                          media_type == 'gal' ~ 'Galactose',
                          media_type == 'spent' ~ 'Filtered E spent media')
    )


media_levels <- c('Galactose', 'Acetate', 'Filtered E spent media', 'Coculture with E')

A <- OD_lab2 %>%
    filter(media_type != "blank") %>% 
    filter(mutation != 'blank') %>% 
    ggplot(aes(x = hour, y = OD , group = well, color = mutation))+
    geom_line()+
    #ggtitle(label = "Growth curves of S WT and S aceA mutant")+
    labs(y = 'Optical density (600 nm)')+
    scale_y_log10()+
    theme_bw(12)+
    labs(color = 'Genotype', x = 'Time (hours)') +
    facet_wrap(~factor(media_type, levels = c('Galactose', 'Acetate', 'Filtered E spent media', 'Coculture with E'))) +
    scale_color_manual(values = c('S WT' = '#F35E5A', 'S∆*aceA*' = '#242F90')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          legend.text = ggtext::element_markdown(color="black"),
          legend.position = c(.63, .99),  
          legend.justification = c(1, 1), 
          axis.text = element_text(color="black"), 
          strip.text = element_text(color="black"),
          axis.title.y = ggtext::element_markdown())

saveRDS(object = A, file = 'rds_plots/aceA_growth_curves.rdata')

# ggsave("plots/SaceA_growth_curves.png",
#        dpi = 300, width = 5, height = 2.5)

## plotting log linear growth rate as a faceted  on mutation and media type

loglin_gr_meta <- loglin_gr %>% left_join(meta, by = "well") %>% 
    separate(values, into = c("mutation", "media_type"), remove = F) %>% 
    mutate(
        mutation = case_when(mutation == 'WT' ~ 'S WT',
                             mutation == 'aceA' ~ 'S∆*aceA*',
                             mutation == 'blank' ~ 'blank'),
        media_type = case_when(media_type == 'ace' ~ 'Acetate',
                               media_type == 'cocult' ~ 'Coculture with E',
                               media_type == 'gal' ~ 'Galactose',
                               media_type == 'spent' ~ 'Filtered E spent media')
    ) %>% 
    full_join(get_yields(OD)) # add final yield too

# growth rate (µ)
B <- loglin_gr_meta %>% 
    filter(media_type != c("NA")) %>%
    filter(media_type != "blank") %>% 
    filter(mutation != "blank") %>% 
    ggplot(aes(x = mutation, y = r, color = mutation)) +
    geom_jitter(shape = 1, width = .05)+
    stat_summary(fun = mean, shape = '-', size = 3, color = 'black')+
    theme_bw(12)+
    theme(legend.position = 'none')+
    labs(x = "", y = "Growth rate (per hour)", color = 'Genotype')+
    facet_wrap(~ factor(media_type, levels = media_levels)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.text.x = ggtext::element_markdown(color="black"), 
          strip.text = element_text(color="black"),
          axis.title.y = ggtext::element_markdown())


# Max OD
C <- loglin_gr_meta %>% 
    filter(media_type != c("NA")) %>%
    filter(media_type != "blank") %>% 
    filter(mutation != "blank") %>% 
    ggplot(aes(x = mutation, y = yield, color = mutation)) +
    geom_jitter(shape = 1, width = .05,)+
    stat_summary(fun = mean, shape = '-', size = 3, color = 'black')+
    theme_bw(12)+
    theme(legend.position = 'none')+
    labs(x = "", y = "Optical density (600 nm)", color = 'Genotype')+
    facet_wrap(~ factor(media_type, levels = media_levels)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.text.x = ggtext::element_markdown(color="black"), 
          strip.text = element_text(color="black"),
          axis.title.y = ggtext::element_markdown())

# Lag time
D <- loglin_gr_meta %>% 
    filter(media_type != c("NA")) %>%
    filter(media_type != "blank") %>% 
    filter(mutation != "blank") %>% 
    ggplot(aes(x = mutation, y = lag, color = mutation)) +
    geom_jitter(shape = 1, width = .05,)+
    stat_summary(fun = mean, shape = '-', size = 3, color = 'black')+
    theme_bw(12)+
    theme(legend.position = 'none')+
    labs(x = "", y = "Lag (hours)", color = 'Genotype')+
    facet_wrap(~ factor(media_type, levels = media_levels)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.text.x = ggtext::element_markdown(color="black"), 
          strip.text = element_text(color="black"),
          axis.title.y = ggtext::element_markdown())


ggpubr::ggarrange(A,B,C,D,common.legend = T,labels = 'AUTO')

ggsave("plots/SaceA_SWT_panels.png",
        dpi = 300, width = 8, height = 6)

#### statistics ####
# t-tests

# compare S∆aceA to S WT growth rate in different conditions

target <- c( "Galactose", "Filtered spent E media", "Coculture with E") #add in "cocult" later
target_mut <- c("S WT","S∆*aceA*")

aceA_WT_mu_ttest <- loglin_gr_meta %>%
    filter(media_type %in% target) %>%
    filter(mutation %in% target_mut) %>%
    group_by(media_type) %>%
    rstatix::t_test(r ~ mutation) %>%
    rstatix::add_significance("p")

#look at lag 
aceA_WT_lag_ttest <- loglin_gr_meta %>%
    filter(media_type %in% target) %>%
    filter(mutation %in% target_mut) %>%
    group_by(media_type) %>%
    rstatix::t_test(lag ~ mutation) %>%
    rstatix::add_significance("p")

#look at yield

# add acetate to the analysis
target <- c("Acetate", "Galactose", "Filtered spent E media", "Coculture with E") #add in "cocult" later
aceA_WT_yield_ttest <- loglin_gr_meta %>% 
    filter(media_type %in% target) %>%
    filter(mutation %in% target_mut) %>%
    group_by(media_type) %>%
    rstatix::t_test(yield ~ mutation) %>%
    rstatix::add_significance("p")
