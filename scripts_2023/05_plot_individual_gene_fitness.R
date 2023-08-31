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
    unite("united_names", c("gene", "gene_name"), sep = " / ", remove = FALSE) %>% 
    mutate(ecology = 
               case_when(ecology == 'mutualism' ~ 'Mutualism',
                         ecology == 'competition' ~ 'Competition'))


##### function for plotting gene fitness ####

gene_fitness <- function(goi, fitness_type, title_var) {  #goi == gene(s) of interest; fitness_type is the where the fitness data is coming from; title_var == title of plot
    fitness_type %>% 
        filter(str_detect(united_names, goi)) %>%
        ggplot(aes(x = factor(community, levels = c('S','SE','SM','SEM')), y = fitness_normalized, color = ecology ))+
        stat_summary(fun = mean, shape = '-', size = 3, color = 'black')+
        geom_point(shape = 1)+
        theme_bw(12)+
        # ggtitle(title_var)+
        xlab("")+
        ylab("Fitness effect")+
        theme( legend.position = 'none') +
        facet_grid(gene_name ~ factor(ecology, levels = c('Mutualism','Competition'))) +
        geom_hline(yintercept = 0, linetype = 'dashed') +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background =element_rect(fill="white"), 
              strip.text.y = element_text(face = "italic", color = 'black'),
              strip.text.x = element_text(color = 'black'),
              axis.text = element_text(color="black") 
              ) # make gene labels italic
}


#### Make figure with nitrogen and galactose gene fitness scores ####

nitrogen <- gene_fitness('amtB|gltD|gltB|lrp|glnK', fit_all, '') #%>% 

# gal <- gene_fitness('galK|galE|galP', fit_all, '')+ theme(legend.position = "none")
galace <- gene_fitness('aceA|galK', fit_all, '')+ theme(legend.position = "none")

nitro_galace_fig<- ggpubr::ggarrange(
    nitrogen,
    galace,
    ncol = 1,
    nrow = 2, 
    heights = c(2,.9)
) 
saveRDS(nitro_galace_fig, 'rds_plots/nitro_galace_barseq_gene_fit.rdata')

# acetate
ace <- gene_fitness('aceA|aceB', fit_all, '')+ theme(legend.position = "none")

##### Mutualism improved fitnesses -- Vitamins and Amino acids #####


simple_amino_vit <- gene_fitness('ilvA|panC|nadC|thiE|pdxB', fit_all, '') +
    facet_grid(factor(gene_name, levels = c('ilvA','panC','nadC','thiE','pdxB')) ~ factor(ecology, levels = c('Mutualism','Competition')))

saveRDS(simple_amino_vit, 'rds_plots/simple_amino_vit_barseq.rdata')

# ggsave(plot = simple_amino_vit, filename = 'plots/simple_amino_vit.pdf', width = 4, height = 8)

