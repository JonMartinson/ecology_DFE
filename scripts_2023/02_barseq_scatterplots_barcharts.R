
rm(list = ls())
library(tidyverse)
library(ggpmisc)

# load data
load("./data/fitness.Rdat")
load("./data/gene_names.Rdat")
stats_mut = read_csv("./data/mutualism_two_way_anova_stats.csv") # mutualism lm stats output for single genes
stats_comp = read_csv("./data/competition_two_way_anova_stats.csv") # competitive lm stats output for single genes
gene_names$locusId[gene_names$gene_name == "ilvA"] = "STM3905"
gene_names$locusId[gene_names$gene_name == "ilvC"] = "STM3909"
gene_names$locusId[gene_names$gene_name == "ilvD"] = "STM3904"
gene_names$locusId[gene_names$gene_name == "ilvE"] = "STM3903"

gene_names = gene_names |>
  rbind(data.frame(gene_name = c("fadA","fadB"),
                   locusId = c("STM3982", "STM3983")))

fitness = fitness %>% left_join(gene_names)

##

# Wrangling data ----------------------------------------------------------

stats_mut <- stats_mut %>% 
  mutate(locusId = gene) %>%
  left_join(gene_names) %>%
    select(term, padj, locusId, gene_name) %>%
    mutate(community = ifelse(term == "(Intercept)",
                              "S",
                              ifelse(term == "ETRUE",
                                     "SE",
                                     ifelse(term == "MTRUE", "SM", "SEM"))))
stats_comp <-  stats_comp %>% 
    mutate(locusId = gene) %>%
    left_join(gene_names)%>%
    select(term, padj, locusId, gene_name) %>%
    mutate(community = ifelse(term == "(Intercept)",
                              "S",
                              ifelse(term == "ETRUE",
                                     "SE",
                                     ifelse(term == "MTRUE", "SM", "SEM"))))
##
# join and summarize the statistic results with the fitness results

fitness_mutualism <- fitness %>%
    filter(ecology == 'mutualism') %>%
    group_by(locusId, E, M, community, gene_name) %>%
    summarize(fitness_mean = mean(fitness_normalized),
              fitness_sd = sd(fitness_normalized) / sqrt(4)) %>%
    left_join(stats_mut) %>%
                  ungroup() %>%
                  select(community, locusId, gene_name, fitness_mean, fitness_sd, padj)  %>%
                 pivot_wider(names_from = "community", values_from = c("fitness_mean", "fitness_sd", "padj")) %>% 
   ungroup()

fitness_competition <- fitness %>%
    filter(ecology == 'competition') %>%
    group_by(locusId, E, M, community, gene_name) %>%
    summarize(fitness_mean = mean(fitness_normalized),
              fitness_sd = sd(fitness_normalized) / sqrt(4)) %>%
    left_join(stats_comp) %>%
    ungroup() %>%
    select(community, locusId, gene_name, fitness_mean, fitness_sd, padj) %>%
    pivot_wider(names_from = "community", values_from = c("fitness_mean", "fitness_sd", "padj")) %>% 
    ungroup()



# plotting functions to make scatter plots --------------------------------


SE_scatter <- function(fitness_df, xlimits = c(-8, 2), ylimits = c(-8, 2), alpha_level = .5){
    fitness_df %>% 
        ggplot(aes(x = fitness_mean_S,
                   y = fitness_mean_SE,
                   color = padj_SE < 0.05))+
        scale_color_brewer(type = "qual", palette = 6, direction = -1)+
        geom_point(data = .%>%
                       filter(padj_SE >= 0.05), alpha = alpha_level)+
        geom_point(data = .%>%
                       filter(padj_SE < 0.05))+
        geom_abline(slope = 1, intercept = 0)+
        theme_bw(12)+
        geom_hline(yintercept = 0, linetype = 'dashed')+
        geom_vline(xintercept = 0, linetype = 'dashed')+
        lims(x = xlimits, y = ylimits) +
        labs(x = "Fitness effect in monoculture",y = "Fitness effect in coculture with E")+
        theme(legend.position = "none",
              panel.grid = element_blank(),
              axis.text = element_text(color="black"), 
              strip.text = element_text(color = 'black'))  
}

SM_scatter <- function(fitness_df, xlimits = c(-8, 2), ylimits = c(-8, 2), alpha_level = .5){
    fitness_df %>% 
        ggplot(aes(x = fitness_mean_S,
                   y = fitness_mean_SM,
                   color = padj_SM < 0.05))+
        scale_color_brewer(type = "qual", palette = 6, direction = -1)+
        geom_point(data = .%>%
                       filter(padj_SM >= 0.05), alpha = alpha_level)+
        geom_point(data = .%>%
                       filter(padj_SM < 0.05))+
        geom_abline(slope = 1, intercept = 0)+
        theme_bw(12)+
        geom_hline(yintercept = 0, linetype = 'dashed')+
        geom_vline(xintercept = 0, linetype = 'dashed')+
        lims(x = xlimits, y = ylimits) +
        labs(x = "Fitness effect in monoculture",y = "Fitness effect in coculture with M")+
        theme(legend.position = "none",
              panel.grid = element_blank(),
              axis.text = element_text(color="black"), 
              strip.text = element_text(color = 'black'))  
}

# in order to get the correct colors, I had to make a specific function for this 
# there has to be a better way to do this. 
SM_competition_scatter <- function(fitness_df, xlimits = c(-8, 2), ylimits = c(-8, 2), alpha_level = .5){
    fitness_df %>% 
        ggplot(aes(x = fitness_mean_S,
                   y = fitness_mean_SM,
                   color = locusId == 'STM0005'))+ # have to set all of the colors to the blue color since there are no significant results
        scale_color_brewer(type = "qual", palette = 6, direction = -1)+
        geom_point(data = .%>%
                       filter(padj_SM >= 0.05), alpha = alpha_level)+
        geom_point(data = .%>%
                       filter(padj_SM < 0.05))+
        geom_abline(slope = 1, intercept = 0)+
        theme_bw(12)+
        geom_hline(yintercept = 0, linetype = 'dashed')+
        geom_vline(xintercept = 0, linetype = 'dashed')+
        lims(x = xlimits, y = ylimits) +
        labs(x = "Fitness effect in monoculture",y = "Fitness effect in coculture with M")+
        theme(legend.position = "none",
              panel.grid = element_blank(),
              axis.text = element_text(color="black"), 
              strip.text = element_text(color = 'black'))  
}


# Scatterplots -----------------------------------------------------------

A <- SE_scatter(fitness_mutualism)

B <- SM_scatter(fitness_mutualism)

C <- SE_scatter(fitness_competition)

D <- SM_competition_scatter(fitness_competition) 



multi_panel <- ggpubr::ggarrange(A,B,C,D)

saveRDS(object = multi_panel, file = 'rds_plots/scatterplot_cocult_vs_mono.rdata')

# If you want to make zoomed in plots of the genes with positive fitness effects

# A_zoom <- SE_scatter(fitness_mutualism, xlimits = c(-1.5, 1.5), ylimits = c(-1.5, 1.5))
# B_zoom <- SM_scatter(fitness_mutualism, xlimits = c(-1.5, 1.5), ylimits = c(-1.5, 1.5))
# C_zoom <- SE_scatter(fitness_competition, xlimits = c(-1.5, 1.5), ylimits = c(-1.5, 1.5))
# D_zoom <- SM_competition_scatter(fitness_competition, xlimits = c(-1.5, 1.5), ylimits = c(-1.5, 1.5))
# multi_panel_zoom <- ggpubr::ggarrange(A_zoom,B_zoom,C_zoom,D_zoom)
# saveRDS(object = multi_panel_zoom, file = 'rds_plots/zoomed_scatterplot_cocult_vs_mono.rdata')

# library(grid) -- if you want to make inset plots with the zoomed in plots

# A_with_inset <- A +
#     annotation_custom(grob = ggplotGrob(A_zoom),  xmin = -8, xmax = -6 , ymin = -2.5, ymax = 2)
# B_with_inset <- B +
#     annotation_custom(grob = ggplotGrob(B_zoom),  xmin = -8, xmax = -6 , ymin = -4.5, ymax = 2.5)
# C_with_inset <- C +
#     annotation_custom(grob = ggplotGrob(C_zoom),  xmin = -8, xmax = -6 , ymin = -4.5, ymax = 2.5)
# D_with_inset <- D +
#     annotation_custom(grob = ggplotGrob(D_zoom),  xmin = -8, xmax = -6 , ymin = -4.5, ymax = 2.5)

# barcharts ------------------------------------------------------

stats_mut_bar = read_csv("./data/mutualism_two_way_anova_stats.csv") # mutualism lm stats output for single genes
stats_comp_bar = read_csv("./data/competition_two_way_anova_stats.csv") # competitive lm stats output for single genes

stats_mut_bar = stats_mut_bar %>%  
    mutate(locusId = gene) %>%
    left_join(gene_names)

stats_comp_bar = stats_comp_bar %>% 
    mutate(locusId = gene) %>%
    left_join(gene_names)

# make data frames with the linear model results and the fitness scores -- also include a community composition variable

mut_fitness <- fitness %>%
    filter(ecology == "mutualism") %>%
    group_by(locusId, E, M, community, gene_name) %>%
    summarize(fitness_mean = mean(fitness_normalized ),
              fitness_sd = sd(fitness_normalized ) / sqrt(4)) %>%
    left_join(stats_mut_bar %>%
                  select(term, padj, locusId, gene_name, Estimate) %>%
                  mutate(community = ifelse(term == "(Intercept)",
                                            "S",
                                            ifelse(term == "ETRUE",
                                                   "SE",
                                                   ifelse(term == "MTRUE", "SM", "SEM"))))) %>%
    ungroup() %>%
    select(community, locusId, gene_name, fitness_mean, fitness_sd, Estimate, padj)

# saveRDS(mut_fitness, 'data/mut_fitness.rdata')

comp_fitness <- fitness %>%
    filter(ecology == "competition") %>%
    group_by(locusId, E, M, community, gene_name) %>%
    summarize(fitness_mean = mean(fitness_normalized ),
              fitness_sd = sd(fitness_normalized ) / sqrt(4)) %>%
    left_join(stats_comp_bar %>%
                  select(term, padj, locusId, gene_name, Estimate) %>%
                  mutate(community = ifelse(term == "(Intercept)",
                                            "S",
                                            ifelse(term == "ETRUE",
                                                   "SE",
                                                   ifelse(term == "MTRUE", "SM", "SEM"))))) %>%
    ungroup() %>%
    select(community, locusId, gene_name, fitness_mean, fitness_sd, Estimate, padj)



##### Count number of significant genes based on initial monoculture fitness, direction of fitness effect, and overlap #####


# identify the sign of the fitness for each gene of the monocultures (mutualism and competition) 

mutualism_monoculture_above_zero_df <- mut_fitness %>% 
    group_by(locusId) %>% 
    mutate(monoculture_above_zero = case_when(community == 'S' & fitness_mean > 0 ~ 'Positive monoculture fitness',
                                              community == 'S' & fitness_mean < 0 ~ 'Negative monoculture fitness')) %>% 
    ungroup() %>% 
    filter(community == 'S') %>% 
    select(locusId, monoculture_above_zero)

competition_monoculture_above_zero_df <- comp_fitness %>% 
    group_by(locusId) %>% 
    mutate(monoculture_above_zero = case_when(community == 'S' & fitness_mean > 0 ~ 'Positive monoculture fitness',
                                              community == 'S' & fitness_mean < 0 ~ 'Negative monoculture fitness')) %>% 
    ungroup() %>% 
    filter(community == 'S') %>% 
    select(locusId, monoculture_above_zero)

#######

# Identify overlap in positive and negative genes in mutualism and competition when considering the direction of effect 
# and the initial fitness of the monoculture for each gene

# mutualism
overlap_mut_full <- mut_fitness %>% 
    filter(padj < 0.05) %>% 
    mutate(direction = 
               case_when(Estimate > 0 ~ 'Positive', 
                         Estimate < 0 ~ 'Negative')) %>% 
    filter(community %in% c('SE', 'SM')) %>% 
    group_by(locusId) %>% 
    left_join(mutualism_monoculture_above_zero_df, by = 'locusId') %>% 
    ungroup() %>% 
    # group_by(direction) %>% 
    pivot_wider(names_from = community, values_from = locusId) %>% select(direction, monoculture_above_zero, SE, SM) 

# Identify the numbe of genes that have overlap when considering the monoculture fitness (negative or positive) and the direction of effect
overlap_mut_filtered <- overlap_mut_full %>%
    group_by(direction, monoculture_above_zero) %>% # group by direction of fitness effect (estimate) and whether the monoculture is above zero
    filter(SE %in% SM & SM %in% SE) %>% # filter the locusIds that are the same in both SM and SE
    dplyr::summarise(number_significant = n()/2) %>% # divide by 2 because gene counts are duplicated in filtering process
    mutate(community = 'Overlap')


# competition --> there is NO overlap, therefore, generated a dummy data frame filled with zeros
comp_fitness %>% 
    filter(padj < 0.05) %>% 
    mutate(direction = 
               case_when(Estimate > 0 ~ 'Positive', 
                         Estimate < 0 ~ 'Negative')) %>% 
    filter(community %in% c('SE', 'SM')) %>% 
    group_by(locusId) %>% 
    left_join(competition_monoculture_above_zero_df, by = 'locusId') %>% 
    ungroup() %>% 
    group_by(community,  direction, monoculture_above_zero) %>% 
    summarise(count = n()) # there is NO possible overlap because has no significant genes

# generate a df to use for the overlap in the figure -- they are all 0
overlap_comp_filtered <- tibble(
    community = c('Overlap', 'Overlap','Overlap','Overlap'), 
    direction = c('Negative', 'Positive', 'Negative', 'Positive'),
    monoculture_above_zero = c('Negative monoculture fitness',
                               'Negative monoculture fitness',
                               'Positive monoculture fitness',
                               'Positive monoculture fitness'),
    number_significant = c(0,0,0,0)
)

###

mutualism_barchart <- mut_fitness %>% 
    filter(padj < 0.05) %>% 
    mutate(direction = 
               case_when(Estimate > 0 ~ 'Positive', 
                         Estimate < 0 ~ 'Negative')) %>% 
    group_by(locusId) %>% 
    left_join(mutualism_monoculture_above_zero_df, by = 'locusId') %>% 
    group_by(community, direction, monoculture_above_zero) %>% 
    summarise(number_significant = n(), .groups = 'drop'
    ) %>% 
    filter(community != 'S') %>% 
    filter(community != 'SEM') %>% 
    mutate(community = case_when(community == 'SE' ~ 'E',
                                 community == 'SM' ~ 'M')) %>% 
    bind_rows(overlap_mut_filtered) %>% 
    complete(community = c("E", "M", "Overlap"), # fill in missing values with 0 for plotting purposes
             direction = unique(direction),
             monoculture_above_zero = unique(monoculture_above_zero), 
             fill = list(number_significant = 0)) %>%
    ggplot(aes(x = community, y = number_significant, fill = direction))+
    geom_col(position = 'dodge', size = .5, color = 'black') +
    geom_text(aes(label = number_significant), vjust = -0.5, position = position_dodge(1))+
    scale_fill_manual(values = c('Negative' = '#7C7568', 'Positive' = 'white'))+ # negative is grey, positive is white
    facet_grid(~monoculture_above_zero) +
    labs(x = '', y = 'Genes with significant effect', fill = 'Direction of effect') +
    # ggtitle('Mutualistic interactions') +
    ylim(0,85) +
    theme_bw(12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(strip.background =element_rect(fill="white"),
          axis.text = element_text(color="black"), 
          strip.text = element_text(color = 'black'))


competition_barplot <- comp_fitness %>% 
    filter(padj < 0.05) %>% 
    mutate(direction = 
               case_when(Estimate > 0 ~ 'Positive', 
                         Estimate < 0 ~ 'Negative')) %>% 
    group_by(locusId) %>% 
    left_join(competition_monoculture_above_zero_df, by = 'locusId') %>% 
    group_by(community, direction, monoculture_above_zero) %>% 
    summarise(number_significant = n(), .groups = 'drop'
    ) %>% 
    filter(community != 'S') %>% 
    filter(community != 'SEM') %>% 
    mutate(community = case_when(community == 'SE' ~ 'E',
                                 community == 'SM' ~ 'M')) %>% 
    bind_rows(overlap_comp_filtered) %>% 
    complete(community = c("E", "M", 'Overlap'), # fiill in missing values with 0 for plotting purposes
             direction = unique(direction),
             monoculture_above_zero = unique(monoculture_above_zero), 
             fill = list(number_significant = 0)) %>%
    ggplot(aes(x = community, y = number_significant, fill = direction))+
    geom_col(position = 'dodge', size = .5, color = 'black') +
    geom_text(aes(label = number_significant), vjust = -0.5, position = position_dodge(1))+
    scale_fill_manual(values = c('Negative' = '#7C7568', 'Positive' = 'white'))+ # negative is grey, positive is white
    facet_grid(~monoculture_above_zero) +
    labs(x = '', y = 'Genes with significant effect', fill = 'Direction of effect') +
    ylim(0,85) +
    theme_bw(12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(strip.background =element_rect(fill="white"),
          axis.text = element_text(color="black"), 
          strip.text = element_text(color = 'black'))


fig_2_bar <- ggpubr::ggarrange(mutualism_barchart, competition_barplot, nrow = 2, legend = 'bottom',common.legend = TRUE)
saveRDS(object = fig_2_bar, file = 'rds_plots/fig_2_bar.rdata')
##

# binomial test results ---------------------------------------------------

# make parsible data frames for the binomial tests


mut_fitness_count_summary<- mut_fitness %>% 
    filter(padj < 0.05) %>% 
    mutate(direction = 
               case_when(Estimate > 0 ~ 'Positive', 
                         Estimate < 0 ~ 'Negative')) %>% 
    group_by(locusId) %>% 
    left_join(mutualism_monoculture_above_zero_df, by = 'locusId') %>% 
    group_by(community, direction, monoculture_above_zero) %>% 
    summarise(number_significant = n(), .groups = 'drop'
    ) %>% 
    filter(community != 'S') %>% 
    filter(community != 'SEM') %>% 
    mutate(community = case_when(community == 'SE' ~ 'E',
                                 community == 'SM' ~ 'M')) %>% 
    bind_rows(overlap_mut_filtered) %>% 
    complete(community = c("E", "M", "Overlap"), # fill in missing values with 0 for plotting purposes
             direction = unique(direction),
             monoculture_above_zero = unique(monoculture_above_zero), 
             fill = list(number_significant = 0)) %>% 
    ungroup()

comp_fitness_count_summary<- comp_fitness %>% 
    filter(padj < 0.05) %>% 
    mutate(direction = 
               case_when(Estimate > 0 ~ 'Positive', 
                         Estimate < 0 ~ 'Negative')) %>% 
    group_by(locusId) %>% 
    left_join(competition_monoculture_above_zero_df, by = 'locusId') %>% 
    group_by(community, direction, monoculture_above_zero) %>% 
    summarise(number_significant = n(), .groups = 'drop'
    ) %>% 
    filter(community != 'S') %>% 
    filter(community != 'SEM') %>% 
    mutate(community = case_when(community == 'SE' ~ 'E',
                                 community == 'SM' ~ 'M')) %>% 
    bind_rows(overlap_mut_filtered) %>% 
    complete(community = c("E", "M", "Overlap"), # fill in missing values with 0 for plotting purposes
             direction = unique(direction),
             monoculture_above_zero = unique(monoculture_above_zero), 
             fill = list(number_significant = 0)) %>% 
    ungroup()



# function to run binomial test on data frames
fitness_binomial_test <- function(df, mono_direction, species){
    df %>%
        filter(monoculture_above_zero == str_c(mono_direction, ' monoculture fitness')) %>%
        filter(community == species) %>%
        summarise(
            sig_negative_genes = number_significant[1],
            total_n = sum(number_significant),
            p_value = binom.test(x = number_significant[1], n = sum(number_significant), p = .5)$p.value
        ) %>% glimpse()
}

# binomial tests for mutualistic conditions
fitness_binomial_test(df = mut_fitness_count_summary, mono_direction = 'Negative', species = 'E')
# $ sig_negative_genes <dbl> 16
# $ total_n            <dbl> 73
# $ p_value            <dbl> 1.525658e-06

fitness_binomial_test(df = mut_fitness_count_summary, mono_direction = 'Positive', species = 'E')
# $ sig_negative_genes <dbl> 12
# $ total_n            <dbl> 12
# $ p_value            <dbl> 0.0004882813


fitness_binomial_test(df = mut_fitness_count_summary, mono_direction = 'Negative', species = 'M')
# $ sig_negative_genes <dbl> 11
# $ total_n            <dbl> 89
# $ p_value            <dbl> 1.367383e-13

fitness_binomial_test(df = mut_fitness_count_summary, mono_direction = 'Positive', species = 'M')
# $ sig_negative_genes <dbl> 6
# $ total_n            <dbl> 7
# $ p_value            <dbl> 0.125


# binomial tests for competitive conditions
fitness_binomial_test(df = comp_fitness_count_summary, mono_direction = 'Negative', species = 'E')
# $ sig_negative_genes <dbl> 22
# $ total_n            <dbl> 46
# $ p_value            <dbl> 0.8829959

fitness_binomial_test(df = comp_fitness_count_summary, mono_direction = 'Positive', species = 'E')
# $ sig_negative_genes <dbl> 28
# $ total_n            <dbl> 29
# $ p_value            <dbl> 1.117587e-07

# this comparison does not work since there are no significant genes in competition with M
# fitness_binomial_test(df = comp_fitness_count_summary, mono_direction = 'Negative', species = 'M')

# this comparison does not work since there are no significant genes in competition with M
# fitness_binomial_test(df = comp_fitness_count_summary, mono_direction = 'Positive', species = 'M')


