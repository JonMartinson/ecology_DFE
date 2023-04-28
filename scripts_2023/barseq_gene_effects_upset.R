rm(list = ls())
library(tidyverse)
library(UpSetR)
library(ggupset)

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

mut_fitness <- fitness %>%
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
    select(community, locusId, gene_name, fitness_mean, fitness_sd, Estimate, padj)

comp_fitness <- fitness %>%
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
    select(community, locusId, gene_name, fitness_mean, fitness_sd, Estimate, padj)


####
mut_fitness %>% 
    filter(padj < 0.05) %>% 
    filter(!community %in% c('S')) %>% 
    group_by(community) %>% 
    summarise(count = n())

comp_fitness %>% 
    filter(padj < 0.05) %>% 
    filter(!community %in% c('S')) %>%
    mutate(
        direction = case_when(
            Estimate > 0 ~ 'Positive',
            Estimate < 0 ~ 'Negative'
        )
    ) %>% 
    group_by(community, direction) %>% 
    summarise(count = n())


df <- mut_fitness %>% 
    filter(padj < 0.05) %>% 
    filter(!community %in% c('S', 'SEM'))

SE <- df %>% filter(community == 'SE')
SM <- df %>% filter(community == 'SM')

all_interact <- bind_rows(SE, SM)

all_sig_list <- list( #create a list with all of the genes that have significance 
    E = SE$locusId,
    M = SM$locusId
)

upset(fromList(all_sig_list), order.by = "freq")
####

SE_positive <- df %>% filter(community == 'SE') %>% filter(Estimate > 0)
SM_positive <- df %>% filter(community == 'SM') %>% filter(Estimate > 0)
SE_negative <- df %>% filter(community == 'SE') %>% filter(Estimate < 0)
SM_negative <- df %>% filter(community == 'SM') %>% filter(Estimate < 0)

all_interact <- bind_rows(SE_positive, SM_positive, SE_negative, SM_negative)

all_sig_list <- list( #create a list with all of the genes that have significance 
    `E +` = SE_positive$locusId,
    `M +` = SM_positive$locusId,
    `E -` = SE_negative$locusId,
    `M -` = SM_negative$locusId
)



upset(fromList(all_sig_list), order.by = "freq", text.scale = 1.25)

#### use ggupset for more options ####

# Combine the data into a single data frame
combined_df <- bind_rows(
    SE_positive %>% mutate(species = "E +"),
    SM_positive %>% mutate(species = "M +"),
    SE_negative %>% mutate(species = "E  -"),
    SM_negative %>% mutate(species = "M  -")
)




combined_df %>% 
    group_by(locusId) %>% 
    summarise(effect = list(species)) %>% 
    ggplot(aes(x = effect)) +
    geom_bar(fill = color_bar) +
    scale_x_upset() +
    theme_classic(12) 
    




###

# Combine the data into a single data frame
combined_df <- bind_rows(
    SE_positive %>% mutate(species = "E +"),
    SM_positive %>% mutate(species = "M +"),
    SE_negative %>% mutate(species = "E  -"),
    SM_negative %>% mutate(species = "M  -")
)

# Create a new column 'combination' to represent different combinations of species
combined_df <- combined_df %>%
    group_by(locusId) %>%
    summarise(effect = list(species)) %>%
    mutate(combination = sapply(effect, function(x) paste(sort(unlist(x)), collapse = " & ")))

# Define custom colors for each combination
color_bar <- c(
    "E +" = "grey", "M +" = "grey", "E  -" = "grey", "M  -" = "grey",
    "E + & M  -" = "grey",
    "E  - & M  -" = "red",
    "E  - & M +" = "grey",
    "E + & M +" = "blue"
)

# Create a bar plot with custom colors
combined_df %>%
    ggplot(aes(x = effect, fill = combination)) +
    geom_bar() +
    scale_x_upset() +
    scale_fill_manual(values = color_bar) +
    theme_classic(12) +
    theme(legend.position = 'none')


# make a bar chart based on upset data

simple_significant_direction<- tribble(
    ~ Species, ~ Direction, ~Number_of_Genes, ~Ecology,
    'E', 'Higher', 15, "Mutualism",
    'E', 'Lower', 20,"Mutualism",
    'M', 'Higher', 36,"Mutualism",
    'M', 'Lower', 10,"Mutualism",
    'Overlap', 'Higher', 37,"Mutualism",
    'Overlap', 'Lower', 2, "Mutualism",
    'E', 'Higher', 50, "Competition",
    'E', 'Lower', 25,"Competition",
    'M', 'Higher', 0,"Competition",
    'M', 'Lower', 0,"Competition",
    'Overlap', 'Higher', 0,"Competition",
    'Overlap', 'Lower', 0, "Competition"
)


simple_significant_direction_will_numbers<- tribble(
    ~ Species, ~ Direction, ~Number_of_Genes, ~Ecology,
    'E', 'Higher', 37+15+5, "Mutualism",
    'E', 'Lower', 20 + 6 +2,"Mutualism",
    'M', 'Higher', 37+ 36+7,"Mutualism",
    'M', 'Lower', 10 +5 +2,"Mutualism",
    'Overlap', 'Higher', 37,"Mutualism",
    'Overlap', 'Lower', 2, "Mutualism",
    'E', 'Higher', 50, "Competition",
    'E', 'Lower', 25,"Competition",
    'M', 'Higher', 0,"Competition",
    'M', 'Lower', 0,"Competition",
    'Overlap', 'Higher', 0,"Competition",
    'Overlap', 'Lower', 0, "Competition"
)
simple_significant_direction %>% 
    ggplot(aes(x = Species, y = Number_of_Genes, fill = Direction )) +
    geom_bar(stat = 'identity', color = 'black', position = position_dodge()) +
    facet_wrap(~factor(Ecology, levels = c('Mutualism', 'Competition')))

overlap_sig_effect_bar<- simple_significant_direction_will_numbers %>% 
    ggplot(aes(x = Species, y = Number_of_Genes, fill = factor(Direction, levels = c('Lower','Higher')) )) +
    geom_bar(stat = 'identity', color = 'black', position = position_dodge()) +
    facet_wrap(~factor(Ecology, levels = c('Mutualism', 'Competition')), ncol = 1) +
    theme_bw(12) +
    # guides(fill = guide_legend(title = 'Direction of Effect')) +
    guides(fill = guide_legend(title = '')) +
    labs(x = '', y = 'Number of Genes \n with Significant Change') +
    theme(legend.position = 'bottom') +
    # scale_fill_manual(values = c('Lower' = 'grey', 'Higher' = 'white'))
    # scale_fill_manual(values = c('Lower' = '#7C7568', 'Higher' = '#D88A42'))
    scale_fill_manual(values = c('Lower' = '#7C7568', 'Higher' = 'white'))

saveRDS(overlap_sig_effect_bar, file = 'rds_plots/overlap_sig_effect_bar.rdata')


mut_overlap_bar <- simple_significant_direction_will_numbers %>% 
    filter(Ecology == 'Mutualism') %>% 
    ggplot(aes(x = Species, y = Number_of_Genes, fill = factor(Direction, levels = c('Lower','Higher')) )) +
    geom_bar(stat = 'identity', color = 'black', position = position_dodge()) +
    facet_wrap(~factor(Ecology, levels = c('Mutualism', 'Competition')), ncol = 1) +
    theme_bw(12) +
    # guides(fill = guide_legend(title = 'Direction of Effect')) +
    guides(fill = guide_legend(title = '')) +
    labs(x = '', y = 'Number of Genes \n with Significant Change') +
    theme(legend.position = 'bottom') +
    # scale_fill_manual(values = c('Lower' = 'grey', 'Higher' = 'white'))
    # scale_fill_manual(values = c('Lower' = '#7C7568', 'Higher' = '#D88A42'))
    scale_fill_manual(values = c('Lower' = '#7C7568', 'Higher' = 'white')) +
    theme(legend.position = 'none')


comp_overlap_bar<- simple_significant_direction_will_numbers %>% 
    filter(Ecology == 'Competition') %>% 
    ggplot(aes(x = Species, y = Number_of_Genes, fill = factor(Direction, levels = c('Lower','Higher')) )) +
    geom_bar(stat = 'identity', color = 'black', position = position_dodge()) +
    facet_wrap(~factor(Ecology, levels = c('Mutualism', 'Competition')), ncol = 1) +
    theme_bw(12) +
    # guides(fill = guide_legend(title = 'Direction of Effect')) +
    guides(fill = guide_legend(title = '')) +
    labs(x = '', y = 'Number of Genes \n with Significant Change') +
    theme(legend.position = 'bottom') +
    # scale_fill_manual(values = c('Lower' = 'grey', 'Higher' = 'white'))
    # scale_fill_manual(values = c('Lower' = '#7C7568', 'Higher' = '#D88A42'))
    scale_fill_manual(values = c('Lower' = '#7C7568', 'Higher' = 'white'))+
    theme(legend.position = 'none')

sig_gene_overlap <- ggpubr::ggarrange(mut_overlap_bar,comp_overlap_bar,
                  ncol = 1, common.legend = F) 

saveRDS(sig_gene_overlap, file = 'rds_plots/sig_gene_overlap.rdata')
