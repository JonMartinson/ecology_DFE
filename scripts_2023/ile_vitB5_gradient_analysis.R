library(tidyverse)
rm(list = ls())

source("scripts_2023/tecan_analysis_functions/tecan_funcs.r")

# function to laod Tecan OD values, subtract the blanks, find the max OD, and join the meta data (well position ids)
yield_map_od <- function(map_path, OD_path){
    meta <- plater::read_plate(file = map_path, #convert 96 well plate map into tidy dataset
                               well_ids_column = "well", #map data has well position identities
                               sep = ",")
    meta$well <- gsub("(?<![0-9])0+", "", meta$well, perl = TRUE) # take out padding 0 in well numbers (A01 --> A1)
    
    # read in the tecan data
    OD <- read_tecan_csv(OD_path, sep = ",", measure = "OD")
    
    # subtracting the blanks
    OD <-  OD %>% group_by(well) %>% 
        mutate(OD = OD - min(OD)) %>%
        mutate(OD = ifelse(OD <= 0.5/200, 0.5/200, OD)) %>%
        ungroup()
    
    OD_yield <- get_yields(OD) # get the max OD value from tecan data
    
    inner_join(OD_yield, meta, "well") %>%  # join meta data and od
    # separate meta data into the mutant strain (mutation), media_type, and concentration of substance
    # that was ammended to media
    separate(col = values, into = c("mutation", "media_type","conc"), sep = "_", remove = F) %>% 
        filter(conc != "NA") %>% 
        mutate(conc = as.numeric(conc))
    }

# read in isoleucine gradient data
ile <- yield_map_od(map_path = 'data/ile_gradient_20221210_map.csv', OD_path = 'data/ile_gradient_20221210_full.csv')

# read in vitamin B5 (pantothenate) gradient data
vitB5 <- yield_map_od(map_path = 'data/vit_B5_big_gradient_map_20221206.csv', OD_path = 'data/B5_gradient_tecan_20221206_complete.csv')

# plot max OD

ile %>%
    filter(media_type != "NA") %>% 
    filter(media_type != "blank") %>% 
    ggplot(aes(x = conc, y= yield, color = media_type)) +
    geom_point(color = "black")+
    geom_point(size = .3)+
    stat_summary(fun=mean, geom="line")+
    labs(x = "Isoleucine (mM)", y = "Max OD600") +
    theme_bw(12)

vitB5 %>%
    filter(media_type != "NA") %>% 
    filter(media_type != "blank") %>% 
    ggplot(aes(x = conc, y= yield, color = media_type)) +
    geom_point(color = "black")+
    geom_point(size = .3)+
    stat_summary(fun=mean, geom="line")+
    labs(x = "Vitamin B5 (µM)", y = "Max OD600") +
    theme_bw(12)

# normalize max OD

yield_normal <- function(gradient){
    gradient %>% 
    group_by(media_type) %>% 
        mutate(adj_yield = yield/max(yield)) %>% 
        ungroup()
}

ile_normal <- yield_normal(ile)
vitB5_normal <- yield_normal(vitB5)

A <- ile_normal %>%
    filter(media_type != "NA") %>%
    filter(media_type != "blank") %>%
    ggplot(aes(x = conc, y = adj_yield, color = media_type)) +
    geom_point(color = "black", size = .6) +
    geom_point(size = .03) +
    stat_summary(fun = mean, geom = "line") +
    labs(x = "Isoleucine (mM)", y = "Normalized Yield", color = 'Carbon Source') +
    theme_bw(12) +
    scale_x_sqrt() +
    scale_color_hue(labels = c('gal' = 'Galactose', 'succ' = 'Succinate'))

B <- vitB5_normal %>%
    filter(media_type != "NA") %>%
    filter(media_type != "blank") %>%
    ggplot(aes(x = conc, y = adj_yield, color = media_type)) +
    geom_point(color = "black", size = .6) +
    geom_point(size = .03) +
    stat_summary(fun = mean, geom = "line") +
    labs(x = "Vitamin B5 (µM)", y = "", color = 'Carbon Source')+
    theme_bw(12) +
    scale_x_sqrt() +
    scale_color_hue(labels = c('gal' = 'Galactose', 'succ' = 'Succinate')) 


ggpubr::ggarrange(A, B, common.legend = T,labels = 'AUTO')

ggsave('plots/vitB5_ile_gradient_normalized_yield.png',
       dpi = 300, height = 2.6, width = 6)
