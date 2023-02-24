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
                       nPerm = 10000,
                       minGSSize = 3,
                       maxGSSize = 800,
                       pvalueCutoff = 0.05,
                       verbose = TRUE,
                       organism = "stm", # 'stm' is the kegg abbreviation for Salmonella enterica LT2
                       pAdjustMethod = "BH")
        
        gse_df <- as.data.frame(gse@result) # convert the gse results to a data frame
        
        gse_df <- gse_df %>%
            mutate(
                species_interact = str_c(species,'_',eco) # create a new column that specifies what 
            )
    }

# Perform gse on each term of the model 

gse_E_mut <- gse_fun(data = df_anova_mut, species = 'ETRUE', eco = 'mutualism')
gse_M_mut <- gse_fun(data = df_anova_mut, species = 'MTRUE', eco = 'mutualism')
gse_EM_mut <- gse_fun(data = df_anova_mut, species = 'ETRUE:MTRUE', eco = 'mutualism') # the E:M interaction is not particularly useful for gse 


# combine the outputs from all of the gse analyses

all_mut <- bind_rows(gse_E_mut, gse_M_mut) #

# all_mut <- bind_rows(gse_E_mut, gse_M_mut,gse_EM_mut)

ggplot(all_mut, aes(x = species_interact,y = fct_reorder(Description, enrichmentScore), fill = (enrichmentScore))) +
    geom_tile() +
    scale_fill_gradient2()+
    theme_minimal() +
    theme(panel.background = element_rect(fill = "gray90"))+
    labs(x = "", y = "Pathway description", fill = "Enrichment Score")

ggsave("plots/gse_heatmap_mutualism.png",
       dpi = 300, width = 8, height = 8)

##### competitive #####
gse_E_comp <- gse_fun(data = df_anova_comp, species = 'ETRUE', eco = 'competition')
gse_M_comp <- gse_fun(data = df_anova_comp, species = 'MTRUE', eco = 'competition')
gse_EM_comp <- gse_fun(data = df_anova_comp, species = 'ETRUE:MTRUE', eco = 'competition') # no pathways are enriched


all_comp <- bind_rows(gse_E_comp, gse_M_comp)


ggplot(all_comp, aes(x = species_interact,y = fct_reorder(Description, enrichmentScore), fill = (enrichmentScore))) +
    geom_tile() +
    scale_fill_gradient2()+
    theme_minimal() +
    theme(panel.background = element_rect(fill = "gray90"))+
    labs(x = "", y = "Pathway description", fill = "Enrichment Score")

ggsave("plots/gse_heatmap_competitition.png",
       dpi = 300, width = 8, height = 2.64)
