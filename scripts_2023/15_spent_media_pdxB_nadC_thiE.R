library(tidyverse)

rm(list = ls())


#grab growth curve plotting functions from Chacon's script 

source("scripts_2023/tecan_analysis_functions/tecan_funcs.r")

#read map data from 96 well plate format
meta <- plater::read_plate(file = "data/spent_media_map_20230818.csv", #convert 96 well plate map into tidy dataset
                           well_ids_column = "well",
                           sep = ",")

meta$well <- gsub("(?<![0-9])0+", "", meta$well, perl = TRUE) # take out padding 0 in well numbers (A01 --> A1) 

# read in the tecan data -- note that this is data from a single time point (it was shaken over 
# a ~3 minute period to disperse cells) -- shaken in an orbital shaker

OD <- read_tecan_csv("data/20230819_BEI_mutant_spent_media_OD600.csv", sep = ",", measure = "OD")

# join OD data with metadata
OD_lab <- inner_join(OD, meta, "well")

# get the od yield

od_yield <- get_yields(OD_lab) %>% 
    left_join(meta, by = "well") %>%
    filter(values != 'na') %>% 
    separate(col = values, into = c("media","strain"), sep = "_", remove = F)


od_yield_summary <- od_yield %>% 
    group_by(strain, media) %>% 
    summarise(
        mean_yield = mean(yield),
        sd_yield = sd(yield)
    ) %>% 
    ungroup()

# calculate the mean blank media value to subtract from the yield for the 
# treatments with cells
mean_blank <- od_yield_summary %>% 
    filter(strain == 'blank') %>% 
    summarise(mean_blank = mean(mean_yield)) %>% 
    pull(mean_blank)

od_yield_summary <- od_yield_summary %>% 
    mutate(mean_yield = mean_yield - mean_blank)

### Plotting 

# background removed blank data plotted:

od_yield %>%
    ggplot(aes(x = strain , y = yield , fill = media)) +
    geom_col(position = 'dodge')




# Update strain names and media labeling information
od_yield_summary_2 <- od_yield_summary %>%
    mutate(
        strain_2 = case_when(
            strain == 'S' ~ 'S',
            strain == 'blank' ~ 'blank',
            TRUE ~ str_c('BEI ∆*', strain, '*')
        ),
        media_2 = case_when(
            media == 'Espent' ~ 'Filtered E spent media',
            media == 'Mspent' ~ 'Filtered M spent media',
            media == 'gal' ~ 'Galactose media',
            media == 'lacmetgal' ~ 'Lactose + galactose + methionine media',
            media == 'succmetgal' ~ 'Succinate + galactose + methionine media',
            TRUE ~ 'Other'  # Default value
        )
    )

# Update strain_level with new names
# strain_level_updated <- c('S', 'blank', str_c('BEI ∆*', c('pdxB', 'nadC', 'thiE'), '*'))
# 
# # Update media_level with new names
# media_level_updated <- c(
#     'Filtered E spent media',
#     'Filtered M spent media',
#     'Galactose media',
#     'Lactose + galactose + methionine media',
#     'Succinate + galactose + methionine media',
#     'Other'  # Default value
# )

# because this is being saved as an rdata file, I can't include the level information as it is not contained in
# the ggplot rdata file. I have to hardcode the levels within the ggplot object. 


BEI_spent_media <- od_yield_summary_2 %>%
    filter(strain != 'blank') %>% 
    ggplot(aes(x = factor(strain_2, levels = c('S', 'blank', str_c('BEI ∆*', c('pdxB', 'nadC', 'thiE'), '*'))) , 
               y = mean_yield , 
               fill = factor(media_2, levels = c(
                   'Filtered E spent media',
                   'Filtered M spent media',
                   'Galactose media',
                   'Lactose + galactose + methionine media',
                   'Succinate + galactose + methionine media',
                   'Other'  # Default value
               )))) +
    geom_errorbar(aes(ymin = mean_yield, ymax = mean_yield + sd_yield), position = position_dodge(0.9), width = 0.5)+
    geom_bar(stat = 'identity', color = 'black', position = position_dodge()) +
    theme_bw(12) +
    labs(x = '', y = 'Absorbance (OD~600~)', fill = 'Media type') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.text.x = ggtext::element_markdown(color="black"), 
          axis.text = element_text(color="black"),
          strip.text = element_text(color="black"),
          axis.title.y = ggtext::element_markdown(color = 'black'),
          axis.title.x = ggtext::element_markdown(color = 'black'),
          legend.position = 'bottom') +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
    scale_fill_viridis_d()

saveRDS(BEI_spent_media, file = 'rds_plots/BEI_spent_media.rdata')


