rm(list = ls())

library(tidyverse)
library(clusterProfiler)

df_anova <- read_csv("data/mutualism_two_way_anova_stats.csv")

## gse on E
df_Etrue <-  df_anova %>% 
  filter(term == "ETRUE")

df_trim_E <- df_Etrue %>% 
  select(c(gene, Estimate))

gene_list_E <-na.omit(df_trim_E)

gene_list_E <- sort(gene_list_E$Estimate, decreasing = TRUE)

original_gene_list_E <- df_trim_E$Estimate

names(original_gene_list_E) <- df_trim_E$gene

gene_list_E<-na.omit(original_gene_list_E)

gene_list_E <- sort(gene_list_E, decreasing = TRUE)


# Run gseKEGG() function to get gene set enrichment results
gse_E <- gseKEGG(geneList = gene_list_E,
                 keyType = "kegg",
                 nPerm = 10000,
                 minGSSize = 3,
                 maxGSSize = 800,
                 pvalueCutoff = 0.05,
                 verbose = TRUE,
                 organism = "stm",
                 pAdjustMethod = "BH")


gse_df_E <- as.data.frame(gse_E@result)

gse_df_E <- gse_df_E %>%
  mutate(
    direction = case_when(
      enrichmentScore > 0 ~ 'Fitness Increased',
      enrichmentScore < 0 ~ 'Fitness Decreased'
    ),
    species_interact = 'Species: E'
  )

ggplot(gse_df_E, aes(x = species_interact, y  = Description, fill = (enrichmentScore))) +
  geom_tile() +
  scale_fill_gradient2()+
  theme_minimal() +
  labs(x = "KEGG Pathway ID", y = "Pathway description", fill = "Enrichment Score") +
  ggtitle("Gene Set Enrichment")


###


df_Mtrue <-  df_anova %>% 
  filter(term == "MTRUE")

df_trim_M <- df_Mtrue %>% 
  select(c(gene, Estimate))

gene_list_M <-na.omit(df_trim_M)

gene_list_M <- sort(gene_list_M$Estimate, decreasing = TRUE)

original_gene_list_M <- df_trim_M$Estimate

names(original_gene_list_M) <- df_trim_M$gene

gene_list_M<-na.omit(original_gene_list_M)

gene_list_M <- sort(gene_list_M, decreasing = TRUE)



# Run gseKEGG() function to get gene set enrichment results
gse_M <- gseKEGG(geneList = gene_list_M,
               keyType = "kegg",
               nPerm = 10000,
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               verbose = TRUE,
               organism = "stm",
               pAdjustMethod = "BH")

# Reshape the gse object to long format for plotting
gse_df_M <- as.data.frame(gse_M@result)

gse_df_M <- gse_df_M %>%
  mutate(
    direction = case_when(
      enrichmentScore > 0 ~ 'Fitness Increased',
      enrichmentScore < 0 ~ 'Fitness Decreased'
    ),
    species_interact = 'Species: M'
    )


# Create heatmap with ggplot2
ggplot(gse_df_M, aes(x = direction, y  = Description, fill = (enrichmentScore))) +
  geom_tile() +
  scale_fill_gradient2()+
  theme_minimal() +
  labs(x = "KEGG Pathway ID", y = "Pathway description", fill = "Enrichment Score") +
  ggtitle("Gene Set Enrichment")

####

ggplot(gse_df_M, aes(x = species_interact, y  = Description, fill = (enrichmentScore))) +
  geom_tile() +
  scale_fill_gradient2()+
  theme_minimal() +
  labs(x = "KEGG Pathway ID", y = "Pathway description", fill = "Enrichment Score") +
  ggtitle("Gene Set Enrichment")


##

gse_df_E_M_mutualism <- full_join(gse_df_E,gse_df_M)

ggplot(gse_df_E_M_mutualism, aes(x = species_interact, y = Description, fill = enrichmentScore)) +
  geom_tile() +
  scale_fill_gradient2() +
  labs(x = "", y = "Pathway description", fill = "Enrichment Score") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "gray90"))+
    ggtitle("Mutualistic Interaction")
ggsave("plots/gse_heatmap_mutualism.png",
       dpi = 300, width = 8, height = 8)

#### Competition ####

df_anova <- read_csv("data/competition_two_way_anova_stats.csv")

## gse on E
df_Etrue <-  df_anova %>% 
    filter(term == "ETRUE")

df_trim_E <- df_Etrue %>% 
    select(c(gene, Estimate))

gene_list_E <-na.omit(df_trim_E)

gene_list_E <- sort(gene_list_E$Estimate, decreasing = TRUE)

original_gene_list_E <- df_trim_E$Estimate

names(original_gene_list_E) <- df_trim_E$gene

gene_list_E<-na.omit(original_gene_list_E)

gene_list_E <- sort(gene_list_E, decreasing = TRUE)


# Run gseKEGG() function to get gene set enrichment results
gse_E <- gseKEGG(geneList = gene_list_E,
                 keyType = "kegg",
                 nPerm = 10000,
                 minGSSize = 3,
                 maxGSSize = 800,
                 pvalueCutoff = 0.05,
                 verbose = TRUE,
                 organism = "stm",
                 pAdjustMethod = "BH")


gse_df_E <- as.data.frame(gse_E@result)

gse_df_E <- gse_df_E %>%
    mutate(
        direction = case_when(
            enrichmentScore > 0 ~ 'Fitness Increased',
            enrichmentScore < 0 ~ 'Fitness Decreased'
        ),
        species_interact = 'Species: E'
    )

ggplot(gse_df_E, aes(x = species_interact, y  = Description, fill = (enrichmentScore))) +
    geom_tile() +
    scale_fill_gradient2()+
    theme_minimal() +
    labs(x = "KEGG Pathway ID", y = "Pathway description", fill = "Enrichment Score") +
    ggtitle("Gene Set Enrichment")


###


df_Mtrue <-  df_anova %>% 
    filter(term == "MTRUE")

df_trim_M <- df_Mtrue %>% 
    select(c(gene, Estimate))

gene_list_M <-na.omit(df_trim_M)

gene_list_M <- sort(gene_list_M$Estimate, decreasing = TRUE)

original_gene_list_M <- df_trim_M$Estimate

names(original_gene_list_M) <- df_trim_M$gene

gene_list_M<-na.omit(original_gene_list_M)

gene_list_M <- sort(gene_list_M, decreasing = TRUE)


# Run gseKEGG() function to get gene set enrichment results
gse_M <- gseKEGG(geneList = gene_list_M,
                 keyType = "kegg",
                 nPerm = 10000,
                 minGSSize = 3,
                 maxGSSize = 800,
                 pvalueCutoff = 0.05,
                 verbose = TRUE,
                 organism = "stm",
                 pAdjustMethod = "BH")

# Reshape the gse object to long format for plotting
gse_df_M <- as.data.frame(gse_M@result)

gse_df_M <- gse_df_M %>%
    mutate(
        direction = case_when(
            enrichmentScore > 0 ~ 'Fitness Increased',
            enrichmentScore < 0 ~ 'Fitness Decreased'
        ),
        species_interact = 'Species: M'
    )


# Create heatmap with ggplot2
ggplot(gse_df_M, aes(x = direction, y  = Description, fill = (enrichmentScore))) +
    geom_tile() +
    scale_fill_gradient2()+
    theme_minimal() +
    labs(x = "KEGG Pathway ID", y = "Pathway description", fill = "Enrichment Score") +
    ggtitle("Gene Set Enrichment")

####

ggplot(gse_df_M, aes(x = species_interact, y  = Description, fill = (enrichmentScore))) +
    geom_tile() +
    scale_fill_gradient2()+
    theme_minimal() +
    labs(x = "KEGG Pathway ID", y = "Pathway description", fill = "Enrichment Score") +
    ggtitle("Gene Set Enrichment")


##

gse_df_E_M_competition <- full_join(gse_df_E,gse_df_M)

ggplot(gse_df_E_M_competition, aes(x = species_interact, y = Description, fill = enrichmentScore)) +
    geom_tile() +
    scale_fill_gradient2() +
    labs(x = "", y = "Pathway description", fill = "Enrichment Score") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "gray90"))+
    ggtitle("Competitive Interaction")

ggsave("plots/gse_heatmap_competitition.png",
       dpi = 300, width = 8, height = 8)

###### 

