rm(list = ls())
library(tidyverse)

# load growth rate data 

growth_rates <- readRDS('data/baranyi_lawn_growth_rates.rdata')

# load fitness data

fitness = read_csv("data/S0_barseq_fitnesses_20220722.csv")

mutualism <- c("T2 S mut", "T3 SE mut", "T5 SM mut", "T5 SEM mut") # early mutualism samples 
competition = c("T6 S comp", "T6 SE comp", "T6 SM comp", "T6 SEM comp")
all_trts = c(mutualism, competition)

fitness = fitness %>%
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


growth_rate_summarised <- growth_rates %>% 
    separate(well, sep = '_', into = c('community', 'replicate')) %>% 
    group_by(community) %>% 
    summarise(mean_growth_rate = mean(r),
              sd_growth_rate = sd(r))


fitness_summarized <- fitness %>% 
    filter(ecology == 'mutualism') %>% 
    group_by(community, rep) %>% 
    summarise(mean_fitness = mean(fitness_normalized)) %>% 
    ungroup() %>% 
    group_by(community) %>% 
    summarise(fitness_mean_of_reps = mean(mean_fitness),
              fitness_sd_of_reps = sd(mean_fitness)) %>% 
    ungroup()

growth_vs_fitness <- left_join(growth_rate_summarised, fitness_summarized, by = 'community') %>% 
    ggplot(aes(x = mean_growth_rate, y = fitness_mean_of_reps, color = factor(community, levels = c('S','SE','SM', 'SEM')))) +
    geom_point(size = 5) +
    geom_errorbarh(aes(xmin = mean_growth_rate - sd_growth_rate, xmax = mean_growth_rate + sd_growth_rate)) +
    geom_errorbar(aes(ymin = fitness_mean_of_reps - fitness_sd_of_reps, ymax = fitness_mean_of_reps + fitness_sd_of_reps), width = .01) +
    labs(color = 'Community',  
         x = 'Mean growth rate', y = 'Mean fitness effect') +
    theme_bw(12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.text.x = ggtext::element_markdown(color="black"), 
          axis.text = element_text(color="black"),
          strip.text = element_text(color="black"),
          axis.title.y = ggtext::element_markdown(color = 'black'),
          axis.title.x = ggtext::element_markdown(color = 'black'),
          legend.position = 'bottom') +
    scale_color_viridis_d()


df_growth_fitness <- left_join(growth_rate_summarised, fitness_summarized, by = 'community') 
cor.test(df_growth_fitness$mean_growth_rate, df_growth_fitness$fitness_mean_of_reps)


saveRDS(growth_vs_fitness, 'rds_plots/growth_vs_fitness.rdata')    




