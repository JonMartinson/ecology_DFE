rm(list = ls())

library(tidyverse)


# load linear model results for mutualism and competition
df_anova_mut <- read_csv("data/mutualism_two_way_anova_stats.csv")
df_anova_comp <- read_csv("data/competition_two_way_anova_stats.csv")

load("./data/gene_names.Rdat")
gene_desc <- read_tsv('data/genes.GC.txt') %>% select(locusId, desc)
gene_names$locusId[gene_names$gene_name == "ilvA"] = "STM3905"
gene_names$locusId[gene_names$gene_name == "ilvC"] = "STM3909"
gene_names$locusId[gene_names$gene_name == "ilvD"] = "STM3904"
gene_names$locusId[gene_names$gene_name == "ilvE"] = "STM3903"

gene_names <- left_join(gene_names, gene_desc, by = 'locusId') %>% 
    unite(col = 'full_description', c('gene_name','desc'), sep =  ' - ', remove = FALSE )

df_anova_comp<- df_anova_comp %>% left_join(gene_names, by = c('gene' = 'locusId'))

df_anova_mut<- df_anova_mut %>% left_join(gene_names, by = c('gene' = 'locusId'))


# function to plot heatmap for genes
heat_plotter <- function(df, species){
    df %>% 
        filter(padj < 0.05) %>% 
        filter(term == species) %>% 
        filter(!is.na(gene_name)) %>%
        mutate(species_interact = case_when(species == 'ETRUE' ~ 'Effect of E', 
                                            species == 'MTRUE' ~ 'Effect of M')) %>%
        ggplot(aes(x = species_interact, y = fct_reorder(full_description, Estimate), fill = Estimate))+
        geom_tile(color = 'black') +
        scale_fill_gradient2(limits = c(-6, 6))+
        theme_minimal()+
        labs(x ='', y = '') +
        theme_bw(12) + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_rect(fill = "white"),
              legend.text = element_text(color="black"),
              axis.text = element_text(color="black"), 
              strip.text = element_text(color="black"))
}

E_mut_effect <- heat_plotter(df_anova_mut, 'ETRUE') + theme(legend.position = 'none')
M_mut_effect <- heat_plotter(df_anova_mut, 'MTRUE') + theme(legend.position = 'none')

E_comp_effect <- heat_plotter(df_anova_comp, 'ETRUE')
M_comp_effect <- heat_plotter(df_anova_comp, 'MTRUE')

mutualism_gene_heatmap<- ggpubr::ggarrange(E_mut_effect, M_mut_effect) %>% 
    ggpubr::annotate_figure( top = ggpubr::text_grob("Mutualistic", 
                                                     color = "black", face = "bold", size = 14))

saveRDS(mutualism_gene_heatmap, file = 'rds_plots/mutualism_gene_heatmap.rdata')

competition_gene_heatmap<- ggpubr::ggarrange(E_comp_effect,common.legend = F) %>% 
    ggpubr::annotate_figure( top = ggpubr::text_grob("Competitive", 
                                                 color = "black", face = "bold", size = 14))

saveRDS(competition_gene_heatmap, file = 'rds_plots/competition_gene_heatmap.rdata')

