rm(list = ls())

library(tidyverse)
library(clusterProfiler)

# load linear model results for mutualism and competition
df_anova_mut <- read_csv("data/mutualism_two_way_anova_stats.csv")
df_anova_comp <- read_csv("data/competition_two_way_anova_stats.csv")

load("./data/gene_names.Rdat")
gene_names$locusId[gene_names$gene_name == "ilvA"] = "STM3905"
gene_names$locusId[gene_names$gene_name == "ilvC"] = "STM3909"
gene_names$locusId[gene_names$gene_name == "ilvD"] = "STM3904"
gene_names$locusId[gene_names$gene_name == "ilvE"] = "STM3903"

df_anova_comp<- df_anova_comp %>% left_join(gene_names, by = c('gene' = 'locusId'))

df_anova_mut<- df_anova_mut %>% left_join(gene_names, by = c('gene' = 'locusId'))



Etrue <- df_anova_mut %>% 
    filter(padj < 0.05) %>% 
    filter(term == 'ETRUE') %>% 
    filter(!is.na(gene_name)) %>%
    mutate(species_interact = 'ETRUE') %>%
    ggplot(aes(x = species_interact, y = fct_reorder(gene_name, Estimate), fill = Estimate))+
    geom_tile() +
    scale_fill_gradient2()+
    theme_minimal()+
    labs(x ='', y = '') +
    theme(legend.position = 'none')

Mtrue <- df_anova_mut %>% 
    filter(padj < 0.05) %>% 
    filter(term == 'MTRUE') %>% 
    filter(!is.na(gene_name)) %>%
    mutate(species_interact = 'MTRUE') %>%
    ggplot(aes(x = species_interact, y = fct_reorder(gene_name, Estimate), fill = Estimate))+
    geom_tile() +
    scale_fill_gradient2()+
    theme_minimal()+
    labs(x ='', y = '') +
    theme(legend.position = 'none')

mutualism_gene_heatmap<- ggpubr::ggarrange(Etrue, Mtrue,common.legend = F) %>% 
    ggpubr::annotate_figure( top = ggpubr::text_grob("Mutualistic", 
                                                     color = "black", face = "bold", size = 14))

saveRDS(mutualism_gene_heatmap, file = 'rds_plots/mutualism_gene_heatmap.rdata')

Etrue_comp <- df_anova_comp %>% 
    filter(padj < 0.05) %>% 
    filter(term == 'ETRUE') %>% 
    filter(!is.na(gene_name)) %>%
    mutate(species_interact = 'ETRUE') %>%
    ggplot(aes(x = species_interact, y = fct_reorder(gene_name, Estimate), fill = Estimate))+
    geom_tile() +
    scale_fill_gradient2()+
    theme_minimal()+
    labs(x ='', y = '')+
    theme(legend.position = 'none')

Mtrue_comp <- df_anova_comp %>% 
    filter(padj < 0.05) %>% 
    filter(term == 'MTRUE') %>% 
    filter(!is.na(gene_name)) %>%
    mutate(species_interact = 'MTRUE') %>%
    ggplot(aes(x = species_interact, y = fct_reorder(gene_name, Estimate), fill = Estimate))+
    geom_tile() +
    scale_fill_gradient2()+
    theme_minimal()+
    labs(x ='', y = '')+
    theme(legend.position = 'none')

competition_gene_heatmap<- ggpubr::ggarrange(Etrue_comp,common.legend = F) %>% 
    ggpubr::annotate_figure( top = ggpubr::text_grob("Competitive", 
                                                 color = "black", face = "bold", size = 14))

saveRDS(competition_gene_heatmap, file = 'rds_plots/competition_gene_heatmap.rdata')

Emut <- df_anova_mut %>% 
    filter(padj < 0.05) %>% 
    filter(term == 'ETRUE') %>% 
    mutate(species_interact = 'ETRUE')

Mmut <- df_anova_mut %>% 
    filter(padj < 0.05) %>% 
    filter(term == 'MTRUE') %>% 
    mutate(species_interact = 'MTRUE')

allmut <- bind_rows(Emut, Mmut)

allmut %>% 
    filter(!is.na(gene_name)) %>%
    ggplot(aes(x = species_interact, y = fct_reorder(gene_name, Estimate), fill = Estimate))+
    geom_tile() +
    scale_fill_gradient2()+
    theme_minimal()+
    coord_flip()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = '', y = '') +
    facet_wrap(~term)
