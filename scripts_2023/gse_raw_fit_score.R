rm(list = ls())

library(tidyverse)
library(clusterProfiler)

load("./data/fitness.Rdat")

###
### function that performs gene set enrichment on the raw mean fitness scores for each community and ecology 
gse_fun <-
    function(species, eco) {
        # species == 'S', 'SE', 'SEM', 'SM'; eco == 'mutualism' or 'competition'
        fit_sort <- fitness %>%
            group_by(locusId, community, ecology) %>%
            summarise(fit_norm = mean(fitness_normalized)) %>%
            ungroup() %>%
            filter(community == species) %>%
            filter(ecology == eco) %>%
            select(locusId, fit_norm) %>%
            arrange(desc(fit_norm))
        
        gene_list <- fit_sort$fit_norm
        names(gene_list) <- fit_sort$locusId
        gse <- gseKEGG(geneList = gene_list,
                       keyType = "kegg",
                       nPerm = 10000,
                       minGSSize = 3,
                       maxGSSize = 800,
                       pvalueCutoff = 0.05,
                       verbose = TRUE,
                       organism = "stm",
                       pAdjustMethod = "BH")
        
        gse_df <- as.data.frame(gse@result)
        
        gse_df <- gse_df %>%
            mutate(
                species_interact = str_c(species,'_',eco)
            )
        
    }

test_S_mut <- gse_fun(species = 'S', eco = 'mutualism')
test_SE_mut <- gse_fun(species = 'SE', eco = 'mutualism')
test_SM_mut <- gse_fun(species = 'SM', eco = 'mutualism')
test_SEM_mut <- gse_fun(species = 'SEM', eco = 'mutualism')


all_mut_interact <- bind_rows(test_S_mut,test_SE_mut,test_SM_mut,test_SEM_mut)

ggplot(
    all_mut_interact,
    aes(x = species_interact, y = fct_reorder(Description, enrichmentScore), fill = enrichmentScore)) +
    geom_tile() +
    scale_fill_gradient2() +
    labs(x = "", y = "Pathway description", fill = "Enrichment Score") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "gray90")) +
    ggtitle("Mutualistic Interaction")


#### competition ####

test_S_comp <- gse_fun(species = 'S', eco = 'competition')
test_SE_comp <- gse_fun(species = 'SE', eco = 'competition')
test_SM_comp <- gse_fun(species = 'SM', eco = 'competition')
test_SEM_comp <- gse_fun(species = 'SEM', eco = 'competition')


all_comp_interact <- bind_rows(test_S_comp,test_SE_comp,test_SM_comp,test_SEM_comp)

ggplot(
    all_comp_interact,
    aes(x = species_interact, y = fct_reorder(Description, enrichmentScore), fill = enrichmentScore)) +
    geom_tile() +
    scale_fill_gradient2() +
    labs(x = "", y = "Pathway description", fill = "Enrichment Score") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "gray90")) +
    ggtitle("Competitive Interaction")
