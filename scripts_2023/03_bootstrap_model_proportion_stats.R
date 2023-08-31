rm(list = ls())
library(tidyverse)
library(ggpmisc)
library(tidymodels)
library(rstatix)

load("./data/fitness.Rdat")
load("./data/gene_names.Rdat")
stats_mut = read_csv("./data/mutualism_two_way_anova_stats.csv")
stats_comp = read_csv("./data/competition_two_way_anova_stats.csv")
gene_names$locusId[gene_names$gene_name == "ilvA"] = "STM3905"
gene_names$locusId[gene_names$gene_name == "ilvC"] = "STM3909"
gene_names$locusId[gene_names$gene_name == "ilvD"] = "STM3904"
gene_names$locusId[gene_names$gene_name == "ilvE"] = "STM3903"

gene_names = gene_names %>% 
    rbind(data.frame(gene_name = c("fadA","fadB"),
                     locusId = c("STM3982", "STM3983")))

fitness = fitness %>% left_join(gene_names)
stats_mut = stats_mut %>% 
    mutate(locusId = gene) %>%
    left_join(gene_names)

stats_comp = stats_comp %>% 
    mutate(locusId = gene) %>%
    left_join(gene_names)

###

# lm functions 
lm_add <- function(split){
    lm(fitness_mean_SEM ~ additive, data = analysis(split)) %>% 
        glance()
}

lm_ave <- function(split){
    lm(fitness_mean_SEM ~ average, data = analysis(split)) %>% 
        glance()
}

lm_strong <- function(split){
    lm(fitness_mean_SEM ~ strongest, data = analysis(split)) %>% 
        glance()
}

### mutualism

mut_models <- fitness %>%
    filter(ecology == "mutualism") %>%
    group_by(locusId, E, M, community, gene_name) %>%
    summarize(fitness_mean = mean(fitness_normalized ),
              fitness_sd = sd(fitness_normalized ) / sqrt(4)) %>%
    left_join(stats_mut %>%
                  select(term, padj, locusId, gene_name, Estimate) %>%
                  mutate(community = ifelse(term == "(Intercept)",
                                            "S",
                                            ifelse(term == "ETRUE",
                                                   "SE",
                                                   ifelse(term == "MTRUE", "SM", "SEM"))))) %>%
    ungroup() %>%
    select(community, locusId, gene_name, fitness_mean, fitness_sd, Estimate, padj) %>%
    pivot_wider(names_from = "community", values_from = c("fitness_mean", "fitness_sd", "padj", "Estimate")) %>%
    mutate(n_SEM = sum(fitness_mean_SEM > fitness_mean_S),
           n_S = sum(fitness_mean_SEM < fitness_mean_S),
           additive = Estimate_S + Estimate_SE + Estimate_SM,
           average = Estimate_S + (Estimate_SE + Estimate_SM)/2,
           strongest = Estimate_S + ifelse(abs(Estimate_SE) > abs(Estimate_SM),
                                           Estimate_SE, Estimate_SM))

#bootstrap mutmodels 


#####
set.seed(123)

mut_boot <- 
    bootstraps(mut_models, times = 2001, apparent = TRUE)





mut_add <- 
    mut_boot %>% 
    mutate(results = map(splits, lm_add))

mut_add_r<- bind_rows(mut_add$results)%>% mutate(type = 'Additive')


mut_ave <- 
    mut_boot %>% 
    mutate(results = map(splits, lm_ave))

mut_ave_r<- bind_rows(mut_ave$results)%>% mutate(type = 'Average')

mut_strong <- 
    mut_boot %>% 
    mutate(results = map(splits, lm_strong))

mut_strong_r<- bind_rows(mut_strong$results) %>% mutate(type = 'Strongest')


mut_all <- bind_rows(mut_strong_r, mut_add_r, mut_ave_r)


mut_all %>% 
    ggplot(aes(r.squared, color = type, fill = type))+
    geom_histogram(bins = 30) +
    facet_wrap(~factor(type, levels = c('Additive', 'Strongest', 'Average')))

mut_box_boot <- mut_all %>% 
    ggplot(aes(x = factor(type, levels = c('Additive', 'Strongest', 'Average')), y = r.squared, color = factor(type, levels = c('Additive', 'Strongest', 'Average')) ))+
    #geom_histogram(bins = 30) +
    #coord_flip()+
    geom_boxplot() +
    theme_bw(12)+
    theme(legend.position = 'none')+
    labs(x = '', y = 'R squared')+
    lims(y = c(.5, 1))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(strip.background =element_rect(fill="white"),
          axis.text = element_text(color="black"), 
          strip.text = element_text(color = 'black'))
#    facet_wrap(~factor(type, levels = c('Additive', 'Strongest', 'Average')))

saveRDS(object = mut_box_boot, file = 'rds_plots/mut_box_bootstrap.rdata')


####### competititive #####
comp_models <- fitness %>%
    filter(ecology == "competition") %>%
    group_by(locusId, E, M, community, gene_name) %>%
    summarize(fitness_mean = mean(fitness_normalized ),
              fitness_sd = sd(fitness_normalized ) / sqrt(4)) %>%
    left_join(stats_comp %>%
                  select(term, padj, locusId, gene_name, Estimate) %>%
                  mutate(community = ifelse(term == "(Intercept)",
                                            "S",
                                            ifelse(term == "ETRUE",
                                                   "SE",
                                                   ifelse(term == "MTRUE", "SM", "SEM"))))) %>%
    ungroup() %>%
    select(community, locusId, gene_name, fitness_mean, fitness_sd, Estimate, padj) %>%
    pivot_wider(names_from = "community", values_from = c("fitness_mean", "fitness_sd", "padj", "Estimate")) %>%
    mutate(n_SEM = sum(fitness_mean_SEM > fitness_mean_S),
           n_S = sum(fitness_mean_SEM < fitness_mean_S),
           additive = Estimate_S + Estimate_SE + Estimate_SM,
           average = Estimate_S + (Estimate_SE + Estimate_SM)/2,
           strongest = Estimate_S + ifelse(abs(Estimate_SE) > abs(Estimate_SM),
                                           Estimate_SE, Estimate_SM))

#bootstrap competitive models 


#####

comp_boot <- 
    bootstraps(comp_models, times = 2001, apparent = TRUE)



comp_add <- 
    comp_boot %>% 
    mutate(results = map(splits, lm_add))

comp_add_r<- bind_rows(comp_add$results)%>% mutate(type = 'Additive')


comp_ave <- 
    comp_boot %>% 
    mutate(results = map(splits, lm_ave))

comp_ave_r<- bind_rows(comp_ave$results)%>% mutate(type = 'Average')

comp_strong <- 
    comp_boot %>% 
    mutate(results = map(splits, lm_strong))

comp_strong_r<- bind_rows(comp_strong$results) %>% mutate(type = 'Strongest')


comp_all <- bind_rows(comp_strong_r, comp_add_r, comp_ave_r)


comp_all %>% 
    ggplot(aes(r.squared, color = type, fill = type))+
    geom_histogram(bins = 30) +
    facet_wrap(~factor(type, levels = c('Additive', 'Strongest', 'Average')))

comp_box_boot <- comp_all %>% 
    ggplot(aes(x = factor(type, levels = c('Additive', 'Strongest', 'Average')), y = r.squared, color = factor(type, levels = c('Additive', 'Strongest', 'Average')) ))+
    #geom_histogram(bins = 30) +
    #coord_flip()+
    geom_boxplot() +
    theme_bw(12)+
    theme(legend.position = 'none')+
    labs(x = '', y = 'R squared') +
    lims(y = c(.5, 1)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(strip.background =element_rect(fill="white"),
          axis.text = element_text(color="black"), 
          strip.text = element_text(color = 'black'))
#    facet_wrap(~factor(type, levels = c('Additive', 'Strongest', 'Average')))

saveRDS(object = comp_box_boot, file = 'rds_plots/comp_box_bootstrap.rdata')

# r2 side by side

ggpubr::ggarrange(mut_box_boot, comp_box_boot)  


bootstrap_models_combined <- bind_rows(Competition = comp_all, Mutualism = mut_all, .id = "id") %>% 
    ggplot(aes(x = factor(type, levels = c('Additive', 'Strongest', 'Average')), y = r.squared, color = factor(type, levels = c('Additive', 'Strongest', 'Average')) ))+
    geom_boxplot() +
    theme_bw(12)+
    theme(legend.position = 'none')+
    # labs(x = '', y = '*R*<sup>2</sup>') +
    labs(x = '', y = '*R*^2^') +
    lims(y = c(.7, 1)) +
    facet_wrap(~factor(id, levels = c('Mutualism', 'Competition'))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(strip.background =element_rect(fill="white")) +
    theme(axis.title.y = ggtext::element_markdown(),
          axis.text = element_text(color="black"), 
          strip.text = element_text(color = 'black'))


saveRDS(object = bootstrap_models_combined, file = 'rds_plots/bootstrap_models_combined.rdata')

ggsave('plots/bootstrap_models.png',
       dpi = 300, height = 2, width = 4.5)

##### Statistics #####

#### Proportion difference tests ####

# Comparing models with proportions#

# Test whether the difference in R-squared values between the compared
# models is greater than zero
# The code computes the p-value as follows:
#
# 1. Count the number of times the differences (i.e., difference in R-squared values) are greater than 0.
# 2. Divide this count by the total number of replicates (found by taking the maximum value of .$replicate) and adding 1.
# 3. Subtract this fraction from 1 to obtain the proportion-based p-value.


mut_all_rep <- mut_all %>% group_by(type) %>%  mutate(replicate = row_number()) %>% ungroup() 
comp_all_rep <- comp_all %>% group_by(type) %>%  mutate(replicate = row_number()) %>% ungroup() 


# Get differences between bootstrap R2s
mut_difference <- mut_all_rep %>% 
    select(replicate,r.squared,type) %>% 
    pivot_wider(names_from = type, values_from = r.squared) %>% 
    mutate(strongest_minus_average = Strongest - Average,
           additive_minus_strongest = Additive - Strongest,
           additive_minus_average = Additive - Average)

comp_difference <- comp_all_rep %>% 
    select(replicate,r.squared,type) %>% 
    pivot_wider(names_from = type, values_from = r.squared) %>% 
    mutate(strongest_minus_average = Strongest - Average,
           additive_minus_strongest = Additive - Strongest,
           additive_minus_average = Additive - Average)



# Perform test in for models in mutualistic and competive models. 

mut_prop_test <- mut_difference %>% 
    select(-c(Additive, Average, Strongest)) %>% 
    pivot_longer(cols = c(strongest_minus_average, additive_minus_strongest, additive_minus_average), 
                 names_to = 'models',
                 values_to = 'differences') %>% 
    group_by(models) %>% 
    summarise(prop_p_val = (1 - ((sum(differences > 0)))) /(max(.$replicate) + 1))

# models                   prop_p_val
# < chr >                         < dbl >
#     1 additive_minus_average     0.000499
#     2 additive_minus_strongest   0.000499
#     3 strongest_minus_average    0.000499

comp_prop_test <- comp_difference %>% 
    select(-c(Additive, Average, Strongest)) %>% 
    pivot_longer(cols = c(strongest_minus_average, additive_minus_strongest, additive_minus_average), 
                 names_to = 'models',
                 values_to = 'differences') %>% 
    group_by(models) %>% 
    summarise(prop_p_val = (1 - ((sum(differences > 0)))) /(max(.$replicate) + 1))

# models                   prop_p_val
# < chr >                         < dbl >
#     1 additive_minus_average     0.000499
#     2 additive_minus_strongest   0.000499
#     3 strongest_minus_average    0.000499

##### Additivity Plot #####

results = read_csv("./data/mutualism_two_way_anova_stats.csv")


## additivity
additivity_bar_plot <- results %>%
    group_by(gene) %>%
    filter(padj[term == "ETRUE"] < 0.05 &
               padj[term == "MTRUE"] < 0.05) %>%
    ungroup() %>%
    filter(term == "ETRUE:MTRUE") %>%
    mutate(interaction = ifelse(padj >= 0.05,
                                "Additive", "Non-additive")) %>%
    mutate(effect = ifelse(
        interaction == "Additive",
        NA,
        ifelse(Estimate > 0, "Positive", "Negative")
    )) %>%
    ggplot(aes(x = interaction, fill = effect)) +
    geom_bar(colour = 'black') +
    labs(x = "Interaction between E and M",
         y = "Genes") +
    theme_bw(12) +
    scale_fill_manual(values = c('Additive' = 'grey',
                                 'Negative' = 'black',
                                 'Positive' = 'white'),
                      breaks = c('Negative', 'Positive'),
                      guide = guide_legend("Change from\nadditivity"))  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(strip.background =element_rect(fill="white")) +
    theme(axis.text = element_text(color="black"), 
          strip.text = element_text(color = 'black'))

saveRDS(object = additivity_bar_plot, file = 'rds_plots/additivity_bar_plot.rdata')

ggsave("./plots/barseq_mutualism_interaction_type_when_both_main_effects_sig.png",
       dpi = 600, width = 5, height = 2.5)
