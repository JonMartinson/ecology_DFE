rm(list = ls())
library(tidyverse)

gene_names <- read_csv("data/uniprot_s_enterica_gene_names.csv") %>%
    mutate(multi_name = `Gene Names`,
           gene_name = `Gene Names (primary)`) %>%
    mutate(stm_name = str_extract_all(multi_name, "STM\\d{4}")) %>% 
    select(gene_name, stm_name)%>% 
    filter(!str_detect(as.character(stm_name), ','))  # remove genes with multiple names

gene_names$stm_name <- unlist(gene_names$stm_name)

load('data/fitness.Rdat')

fit_all <- fitness %>% 
    rename(gene = locusId) %>% 
    left_join(gene_names, by = c('gene' = 'stm_name') ) %>% 
    unite("united_names", c("gene", "gene_name"), sep = " / ", remove = FALSE)


level_order <-   c('S','SE','SM','SEM') #figure ordering

gene_fitness <- function(goi, fitness_type, title_var) {  #goi == gene(s) of interest; fitness_type is the where the fitness data is coming from; title_var == title of plot
    fitness_type %>% 
        filter(str_detect(united_names, goi)) %>%
        ggplot(aes(x = factor(community, levels = level_order), y = fitness_normalized, color = ecology ))+
        stat_summary(fun = mean, shape = '-', size = 3, color = 'black')+
        geom_point(shape = 1)+
        theme_bw(12)+
        ggtitle(title_var)+
        xlab("")+
        ylab("")+
        theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = 'none') +
        facet_grid(gene_name ~ factor(ecology, levels = c('mutualism','competition')))
}

##### Nitrogen stress genes #####

gene_fitness('amtB|gltD|gltB|lrp|glnK', fit_all, 'Nitrogen stress') %>% 
    ggpubr::annotate_figure(left = ggpubr::text_grob("Fitness Normalized", rot = 90))

ggsave("plots/barseq_nitrogen_gene_fitness.png",
       dpi = 300, width = 8, height = 8)


##### acetate and galactose #####

ace <- gene_fitness('aceA|aceB', fit_all, 'Acetate Metabolism')+ theme(legend.position = "none")
gal <- gene_fitness('galK|galE|galP', fit_all, 'Galactose Metabolism')+ theme(legend.position = "none")

ggpubr::ggarrange(
    ace,
    gal,
    labels = c('A', 'B'),
    ncol = 1,
    nrow = 2
) %>%
    ggpubr::annotate_figure(left = ggpubr::text_grob("Fitness Normalized", rot = 90))

ggsave("plots/barseq_ace_gal_fitness.png",
       dpi = 300, width = 8, height = 8)


##### Mutualism improved fitnesses #####

nad <- gene_fitness('nadB|nadC', fit_all, 'NAD biosynthesis')+ theme(legend.position = "none")
thi <- gene_fitness('thiD|thiE', fit_all, 'Vitamin B1 biosynthesis')+ theme(legend.position = "none")
pan <- gene_fitness('panB|panC', fit_all, 'Vitamin B5 biosynthesis')+ theme(legend.position = "none")
pdx <- gene_fitness('STM0091|pdxB', fit_all, 'Vitamin B6 biosynthesis')+ theme(legend.position = "none")
ilv <- gene_fitness('ilvA|ilvC', fit_all, 'Isoleucine Biosynthesis')


ggpubr::ggarrange(
    nad,
    thi,
    pan,
    pdx,
    ilv,
    labels = c('A', 'B', 'C', 'D', 'E'),
    ncol = 2,
    nrow = 3,
    common.legend = F
) %>%
    ggpubr::annotate_figure(left = ggpubr::text_grob("Fitness Normalized", rot = 90))

ggsave("plots/barseq_mutualistic_fitness_improvement.png",
       dpi = 300, width = 8, height = 8)
