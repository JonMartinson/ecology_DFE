rm(list = ls())

library(tidyverse)

source("scripts_2023/tecan_analysis_functions/tecan_funcs.r")

# Note this data is derived from a single timepoint after XX hours of growth 

#read map data from 96 well plate format
meta <- plater::read_plate(file = "data/S_mutant_map_20230217.csv", #convert 96 well plate map into tidy dataset
                           well_ids_column = "well",
                           sep = ",")


meta$well <- gsub("(?<![0-9])0+", "", meta$well, perl = TRUE) # take out padding 0 in well numbers (A01 --> A1) 

# read in the tecan data

OD <- read_tecan_csv("data/salmonella_mutant_ilvA_panC_20230217.csv", sep = ",", measure = "OD")

OD_lab <- inner_join(OD, meta, "well")




od_yield <- get_yields(OD_lab) %>% 
    left_join(meta, by = "well") %>% 
    separate(col = values, into = c("mutation", "media_type"), sep = "_", remove = F) %>% 
    mutate(
        mutation = case_when(mutation == 'SWT' ~ 'S WT',
                             mutation == 'SpanC' ~ 'S∆*panC*',
                             mutation == 'SilvA' ~ 'S∆*ilvA*',
                             mutation == 'blank' ~ 'blank'),
        media_type = case_when(media_type == 'gal' ~ 'Galactose',
                               media_type == 'galile' ~ 'Galactose + Isoleucine',
                               media_type == 'galB5' ~ 'Galactose + Vitamin B5',
                               media_type == 'Espent' ~ 'Filtered E spent media',
                               media_type == 'Mspent' ~ 'Filtered M spent media')
    )

# get blanks OD values for each media_type
media_blanks <- od_yield %>% 
    group_by(mutation, media_type) %>% 
    filter(mutation == 'blank') %>% 
    filter(!(is.na(media_type))) %>% 
    summarise(
        blank = yield) %>% 
    ungroup() %>% 
    select(!mutation)

od_yield <- od_yield %>% 
    left_join(media_blanks, by = 'media_type') %>% # join media_blank df to od_yield
    mutate(yield_subtract = yield - blank, # subtract blank value for each media type from yield
           yield_subtract = case_when(yield_subtract < 0 ~ 0,# for plotting purposes, if the yield is less than OD 0, make the value 0 since there is no growth
                                      yield_subtract >= 0 ~ yield_subtract)) 

yield_summary <- od_yield %>% 
    group_by(media_type, mutation) %>% 
    summarise(
        mean_yield = mean(yield_subtract),
        sd_yield = sd(yield_subtract)
    )
    

yield_summary %>% 
    filter(mutation != "blank") %>% 
    ggplot(aes(x = mutation, y = mean_yield, fill = mutation))+
    geom_bar(stat = "identity", color = "black", position =  position_dodge())+
    geom_errorbar(aes(ymin=mean_yield, ymax = mean_yield + sd_yield), width = .5)+
    labs(
        x = "", y = "Optical density (600 nm)")+
    facet_wrap( ~factor(media_type, levels = 
                            c("Galactose",
                              "Galactose + Isoleucine",
                              "Galactose + Vitamin B5",
                              "Filtered E spent media",
                              "Filtered M spent media"))) +
    theme_bw(12)+
    theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background =element_rect(fill="white"), 
          axis.text.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown(),
          axis.text = element_text(color="black"),
          strip.text = element_text( color = 'black')
          ) 

ggsave('plots/yield_S_mutants_spent_media.png',
       dpi = 300, height = 2, width = 6)

EM_spent_ilvA_panC<- yield_summary %>% 
    filter(mutation != "blank") %>% 
    filter(!media_type %in% c('Galactose + Isoleucine', 'Galactose + Vitamin B5')) %>% 
    ggplot(aes(x = mutation, y = mean_yield, fill = mutation))+
    geom_bar(stat = "identity", color = "black", position =  position_dodge())+
    geom_errorbar(aes(ymin=mean_yield, ymax = mean_yield + sd_yield), width = .5)+
    labs(
        x = "", y = "Optical density (600 nm)")+
    facet_wrap( ~factor(media_type, levels = 
                            c("Galactose",
                              "Filtered E spent media",
                              "Filtered M spent media"))) +
    theme_bw(12)+
    theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background =element_rect(fill="white"), 
          axis.text.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown(),
          axis.text = element_text(color="black"),
          strip.text = element_text( color = 'black')
    ) 

saveRDS(object = EM_spent_ilvA_panC, file = 'rds_plots/EM_spent_ilvA_panC.rdata')

ggsave('plots/yield_S_mutants_spent_media_no_supplementation.png',
       dpi = 300, height = 2, width = 6)


ilvA_panC_supplemented_yield <- yield_summary %>% 
    filter(mutation != "blank") %>% 
    filter(media_type %in% c('Galactose + Isoleucine', 'Galactose + Vitamin B5')) %>% 
    ggplot(aes(x = mutation, y = mean_yield, fill = mutation))+
    geom_bar(stat = "identity", color = "black", position =  position_dodge())+
    geom_errorbar(aes(ymin=mean_yield, ymax = mean_yield + sd_yield), width = .5)+
    labs(
        x = "", y = "Optical density (600 nm)")+
    facet_wrap( ~factor(media_type, levels = 
                            c("Galactose",
                              'Galactose + Isoleucine',
                              'Galactose + Vitamin B5'))) +
    theme_bw(12)+
    theme(legend.position = "none") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background =element_rect(fill="white"), 
          axis.text.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown(),
          axis.text = element_text(color="black"),
          strip.text = element_text( color = 'black')
    ) 

saveRDS(ilvA_panC_supplemented_yield, file = 'rds_plots/ilvA_panC_supplemented_yield.rdata')
