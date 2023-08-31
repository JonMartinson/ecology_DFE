rm(list = ls())
library(tidyverse)

fitness <- read_csv("data/S0_barseq_fitnesses_20220722.csv")

mutualism <- c("T2 S mut", "T3 SE mut", "T5 SM mut", "T5 SEM mut") # early mutualism samples 
competition <- c("T6 S comp", "T6 SE comp", "T6 SM comp", "T6 SEM comp")
all_trts <- c(mutualism, competition)

fitness <- fitness %>%
    filter(treatment %in% all_trts) %>%
    group_by(treatment, locusId, rep, replicate, desc) %>%
    summarize(fitness_normalized = first(fitness_normalized),
              total_raw_counts = sum(raw_count)) %>%
    ungroup() %>%
    mutate(ecology = ifelse(str_detect(treatment, "mut"), "mutualism", "competition")) %>%
    mutate(ecology = factor(ecology, levels = c("mutualism", "competition"))) %>%
    mutate(E = ifelse(str_detect(treatment, "E"), TRUE, FALSE), # add variable whether E or M is present in the treatment
           M = ifelse(str_detect(treatment, "M "), TRUE, FALSE)) %>%
    mutate(community = ifelse(E & M, "SEM",
                              ifelse(E, "SE",
                                     ifelse(M, "SM", "S")))) 

# calculate the median fitness effect  for each community

fitness_median <- fitness %>%
    mutate(community = factor(community,
                              levels = c("S","SE","SM","SEM"))) %>%
    mutate(ecology = 
               case_when(
                   ecology == 'mutualism' ~ 'Mutualism',
                   ecology == 'competition' ~ 'Competition'
               )) %>% 
    group_by(community, ecology, rep) %>%
    summarize(
        median_fitness = median(fitness_normalized),
        ) %>% 
    ungroup()



median_fitness_effect <- fitness_median %>% 
    ggplot(aes(x = community, y = median_fitness, color = community)) +
    stat_summary(fun = mean, shape = "-", size = 3, color = "black") +
    geom_point(shape = 1)+
    facet_wrap(~factor(ecology, levels = c('Mutualism', 'Competition')))+
    theme_bw(12)+
    scale_color_discrete(guide = "none")+
    #geom_text(aes(label = replicate))+
    labs(x = "Community members",
         y = "Median fitness effect") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          legend.text = element_text(color="black", size = 12),
          axis.text = element_text(color="black", size = 12), 
          strip.text = element_text(color="black", size = 12)) +
    geom_hline(yintercept = 0, linetype = 'dashed')


saveRDS(median_fitness_effect, file = 'rds_plots/median_fitness_effect.rdata')
