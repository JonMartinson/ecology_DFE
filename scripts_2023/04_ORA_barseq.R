# Perform overrepresentation analysis on barseq data

rm(list = ls())

library(tidyverse)
library(clusterProfiler)

# load linear model results for mutualism and competition
df_anova_mut <- read_csv("data/mutualism_two_way_anova_stats.csv")
df_anova_comp <- read_csv("data/competition_two_way_anova_stats.csv")

# function that generates gene set enrichment analysis output

gse_fun <-
    function(data, species, eco) { # data == one of the two data frames above; species == 'ETRUE', 'MTRUE', 'ETRUE:MTRUE'; eco == 'mutualism' or 'competition'
        gene_sort <- data %>%
            filter(term == species) %>% 
            select(gene, Estimate) %>%
            arrange(desc(Estimate)) # sort the genes based on Estimate score
        
        gene_list <- gene_sort$Estimate # make a new list with just the Estimates
        names(gene_list) <- gene_sort$gene # add names to the new list of sorted genes
        
        gse <- gseKEGG(geneList = gene_list, # gse analysis
                       keyType = "kegg",
                       # nPerm = 10000,
                       # minGSSize = 3,
                       # maxGSSize = 800,
                       pvalueCutoff = 0.05,
                       verbose = TRUE,
                       organism = "stm", # 'stm' is the kegg abbreviation for Salmonella enterica LT2
                       pAdjustMethod = "BH")
        
        gse_df <- as.data.frame(gse@result) # convert the gse results to a data frame
        
        gse_df <- gse_df %>%
            mutate(
                species_interact = str_c(species,'_',eco) # create a new column that specifies what 
            ) %>% 
            separate(Description, into = c('Description', 'Species'), sep = ' - ')
        
    }

# Perform gse on each term of the model 

gse_E_mut <- gse_fun(data = df_anova_mut, species = 'ETRUE', eco = 'mutualism')
gse_M_mut <- gse_fun(data = df_anova_mut, species = 'MTRUE', eco = 'mutualism')
gse_EM_mut <- gse_fun(data = df_anova_mut, species = 'ETRUE:MTRUE', eco = 'mutualism') # the E:M interaction is not particularly useful for gse 


# combine the outputs from all of the gse analyses

all_mut <- bind_rows(gse_E_mut, gse_M_mut) #

# all_mut <- bind_rows(gse_E_mut, gse_M_mut,gse_EM_mut)




ggplot(all_mut, aes(x = species_interact,y = fct_reorder(Description, enrichmentScore), fill = (enrichmentScore))) +
    geom_tile(color = "black", size = 0.5) +
    scale_fill_gradient2()+
    theme_minimal() +
    # theme(panel.background = element_rect(fill = "gray90"))+
    labs(x = "", y = "Pathway description", fill = "Enrichment Score", color = 'black')  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(strip.background =element_rect(fill="white"),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          axis.text = element_text(color="black"), 
          strip.text = element_text(color = 'black')) 

ggsave("plots/gse_heatmap_mutualism.png",
       dpi = 300, width = 8, height = 8)

##### competitive #####
gse_E_comp <- gse_fun(data = df_anova_comp, species = 'ETRUE', eco = 'competition')
gse_M_comp <- gse_fun(data = df_anova_comp, species = 'MTRUE', eco = 'competition')
gse_EM_comp <- gse_fun(data = df_anova_comp, species = 'ETRUE:MTRUE', eco = 'competition') # no pathways are enriched


all_comp <- bind_rows(gse_E_comp, gse_M_comp)


ggplot(all_comp, aes(x = species_interact,y = fct_reorder(Description, enrichmentScore), fill = (enrichmentScore)))+
    geom_tile(color = "black", size = 0.5) +
    scale_fill_gradient2()+
    theme_minimal() +
    # theme(panel.background = element_rect(fill = "gray90"))+
    labs(x = "", y = "Pathway description", fill = "Enrichment Score", color = 'black')  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(strip.background =element_rect(fill="white"),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          axis.text = element_text(color="black"), 
          strip.text = element_text(color = 'black')) 

ggsave("plots/gse_heatmap_competitition.png",
       dpi = 300, width = 8, height = 2.64)



#### Overrepresentation analysis ####





# denote pathway description that we want to point out in figure

#### ORA with universe

ora_fun <- function(data, species, eco) {
    sig_direction <- data %>% 
        filter(term == species) %>% 
        filter(padj < 0.05) %>% 
        mutate(direction = case_when(Estimate > 0 ~ 'positive',
                                     Estimate < 0 ~ 'negative'))
    positive <- sig_direction %>% 
        filter(direction == 'positive')
    positive <- positive$gene
    positive <- enrichKEGG(positive, organism = 'stm', universe = (unique(data$gene)), keyType = 'kegg')
    
    
    negative <- sig_direction %>% 
        filter(direction == 'negative')
    negative <- negative$gene
    negative <- enrichKEGG(negative, organism = 'stm', keyType = 'kegg')
    
    positive <- positive@result %>% 
        mutate(direction = 'positive')
    negative <- negative@result %>% 
        mutate(direction = 'negative')
    
    alldir <- bind_rows(positive, negative) %>% 
        separate(Description, into = c('Description', 'Species'), sep = ' - ') %>% 
        filter(p.adjust < 0.05) %>% 
        mutate(species_interact = str_c(species, '_', eco))
}

##

ETRUE_ora <- ora_fun(data = df_anova_mut, species = 'ETRUE', eco = 'mutualism')
MTRUE_ora <- ora_fun(data = df_anova_mut, species = 'MTRUE', eco = 'mutualism')

mut_ora <-  bind_rows(ETRUE_ora, MTRUE_ora)

path_cats <- read_csv('data/stm_kegg_pathway_categories.csv')

path_cats<- path_cats %>%
    mutate(ID = substr(pathway, 1, 5)) %>% 
    mutate(ID = str_c('stm',ID)) %>% 
    mutate(ID_rank = rank(category, ties.method =  'first'))



mut_ora <- mut_ora %>% 
    left_join(path_cats,by = 'ID')




# denote pathway description that we want to point out in figure

special <- c('Valine, leucine and isoleucine biosynthesis', 'Nitrogen metabolism', 'Pantothenate and CoA biosynthesis')

mut_ora <- mut_ora %>% 
    mutate(focus = case_when(Description %in% special ~ 'special',
                             TRUE ~ 'not special'))


    

mut_ora_plot <- mut_ora %>%
    mutate(species_interact = case_when(
        species_interact == 'ETRUE_mutualism' ~ 'E effect',
        species_interact == 'MTRUE_mutualism' ~ 'M effect'),
        direction = case_when(direction == 'positive' ~ '+',
                              direction == 'negative' ~ '-')) %>%
    filter(category != 'Global and overview maps') %>%
    rename(Category = category) %>%
    ggplot(aes(x = direction, y = fct_reorder(Description,desc(ID_rank)), color = Category, fill = Category, shape = focus)) +
    geom_point(size = 4, stroke = 1) +
    scale_shape_manual(values = c(19,10)) +
    labs(y = '', x = 'Direction of effect relative to monoculture') +
    facet_wrap(~species_interact) +
    scale_color_brewer(palette = 2, type = 'qual') +
    theme_bw(12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          legend.position = 'bottom',         # Bottom right corner
          legend.justification = c(0, 1),
          axis.text = element_text(color="black", size = 12), 
          legend.text = element_text(color="black", size = 12),
          strip.text = element_text(color="black", size = 12)) +  # Justification
    guides(shape = FALSE) +
    guides(color = guide_legend(nrow = 5, byrow = TRUE)) +
    scale_color_viridis_d() +
    scale_fill_viridis_d()




saveRDS(object =mut_ora_plot, file = 'rds_plots/mutualism_ORA.rdata')

ggsave("plots/ORA_barseq.png",
       dpi = 300, width =8, height = 4)
