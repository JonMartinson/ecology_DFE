rm(list = ls())
library(tidyverse)
library(ggridges)
options(dplyr.summarise.inform = FALSE)

library(tidymodels)

## this contains the final analysis. Some tests and alternatives are listed below
## some of these plots nevertheless will end up in supp

## test for outlier samples using robust PCA is in ./r_scripts/barseq_outlier_tests.r
## test for effect of changing fitness of genes with zero counts in T1 
##     is in ./r_scripts/barseq_effect_of_diff_t1_zerocount_assumptions.r
## old plots like one-sample t-tests and volcano plots are in barseq_dfe_analysis_old.r

fitness = read_csv("data/S0_barseq_fitnesses_20220722.csv")
gene_names = read_csv("data/uniprot_s_enterica_gene_names.csv") %>%
  mutate(multi_name = `Gene Names`,
         gene_name = `Gene Names (primary)`) %>%
  mutate(locusId = sapply(multi_name,
                           FUN = function(x){
                             str_split(x, " ")[[1]][length(str_split(x, " ")[[1]])]
                           })) %>%
  select(gene_name, locusId)

mutualism <- c("T2 S mut", "T3 SE mut", "T5 SM mut", "T5 SEM mut") # early mutualism samples 
competition = c("T6 S comp", "T6 SE comp", "T6 SM comp", "T6 SEM comp")
all_trts = c(mutualism, competition)

fitness = fitness %>%
  filter(treatment %in% all_trts) %>%
  group_by(treatment, locusId, rep, replicate, desc) %>%
  summarize(fitness_normalized = first(fitness_normalized),
            total_raw_counts = sum(raw_count)) %>%
  ungroup() %>%
  mutate(ecology = ifelse(str_detect(treatment, "mut"), "mutualism", "competition")) %>%
  mutate(ecology = factor(ecology, levels = c("mutualism", "competition"))) %>%
  mutate(E = ifelse(str_detect(treatment, "E"), TRUE, FALSE), # add variable whether E or M is present in the treatment
         M = ifelse(str_detect(treatment, "M "), TRUE, FALSE)) %>%
  mutate(community = ifelse(E & M, "SEM",
                            ifelse(E, "SE",
                                   ifelse(M, "SM", "S")))) 

################
#  Basic DFE plots; hists, violins, summary statistics
###############@

# ridge plots with jittered points below#



DFE_histogram <- fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
    mutate(ecology = 
               case_when(
                   ecology == 'mutualism' ~ 'Mutualism',
                   ecology == 'competition' ~ 'Competition'
               )) %>% 
  group_by(community, ecology, locusId) %>%
  summarize(fitness_normalized = mean(fitness_normalized)) %>%
  ggplot(aes(x = fitness_normalized, fill = community))+
  geom_histogram(binwidth = 0.25)+
  theme_bw(12)+
  scale_y_continuous(trans = "log2")+
  scale_fill_discrete(guide = "none")+
  facet_grid(community~factor(ecology, levels = c('Mutualism', 'Competition')))+
  labs(x = "Gene-level fitness effect",
       y = "Genes") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(strip.background =element_rect(fill="white"),
        axis.text = element_text(color="black"),
        strip.text = element_text(color = 'black'),
        axis.title.y = ggtext::element_markdown(color = 'black'))

saveRDS(object = DFE_histogram, file = 'rds_plots/DFE_histogram.rdata')

ggsave("plots/barseq_dfe_log2_histogram.png",
       dpi = 300, width = 5, height = 2.7)



mean_fit_effect_plot <- fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
    mutate(ecology = 
               case_when(
                   ecology == 'mutualism' ~ 'Mutualism',
                   ecology == 'competition' ~ 'Competition'
               )) %>% 
  group_by(community, ecology, replicate) %>%
  summarize(mean_fitness = mean(fitness_normalized)) %>%
  ggplot(aes(x = community, y = mean_fitness, color = community))+
  stat_summary(fun = mean, shape = "-", size = 3, color = "black")+
  geom_point(shape = 1)+
  facet_wrap(~factor(ecology, levels = c('Mutualism', 'Competition')))+
  theme_bw(12)+
  scale_color_discrete(guide = "none")+
  #geom_text(aes(label = replicate))+
  labs(x = "Community members",
       y = "Mean fitness effect") +
  theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(strip.background =element_rect(fill="white"),
          axis.text = element_text(color="black"),
          strip.text = element_text(color = 'black'))
    

saveRDS(object = mean_fit_effect_plot, file = 'rds_plots/mean_fit_effect_plot.rdata')

ggsave("./plots/barseq_dfe_mean_fitness_effect_mean_stat.png",
       dpi = 300, width = 4, height = 2.5)

#statistics--parametric regression

regressions <- fitness %>%
    mutate(community = factor(community,
                              levels = c("S","SE","SM","SEM"))) %>%
    group_by(community, ecology, replicate, E, M, treatment) %>%
    summarize(mean_fitness = mean(fitness_normalized)) %>%
    ungroup() %>% 
    nest(data = c(-ecology)) %>% 
    mutate(
        fit = map(data, ~ lm(mean_fitness ~ E * M, data = .x)),
        tidied = map(fit, tidy), # make clean df with coefficient terms
        glanced = map(fit, glance), # get R2 values and other valuable stuff
        augmented = map(fit, augment)
    )


tidied_model_fits <- regressions %>% # contains model coefficients, p values, etc
    select(ecology, tidied) %>% 
    unnest(tidied)

write_csv(tidied_model_fits, file = 'data/linear_regression_coefficients.csv') 

tidied_model_glance <- regressions %>% # contains R2, F statistic
    select(ecology, glanced) %>% 
    unnest(glanced)

write_csv(tidied_model_glance, file = 'data/linear_regression_glance_metrics.csv') #

##### Stdev in fitness effects #####
sd_fitness <- fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  group_by(community, ecology, replicate) %>%
  summarize(mean_fitness = sd(fitness_normalized)) %>%
  ggplot(aes(x = community, y = mean_fitness, color = community))+
  stat_summary(fun = mean, shape = "-", size = 3, color = "black")+
  geom_point(shape = 1)+
  facet_wrap(~ecology)+
  theme_bw(12)+
  scale_color_discrete(guide = "none")+
  #geom_text(aes(label = replicate))+
  labs(x = "community members",
       y = "Standard deviation of \n fitness effect")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    theme(strip.background =element_rect(fill="white"),
          axis.text = element_text(color="black"),
          strip.text = element_text(color = 'black'))
    
saveRDS(object = sd_fitness, file = 'rds_plots/sd_fitness.rdata')
    
ggsave("./plots/barseq_dfe_sd_fitness.png",
       dpi = 300, width = 4, height = 2.5)

#statistics on standard deviation
fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  group_by(community, ecology, replicate, E, M, treatment) %>%
  summarize(sd_fitness = sd(fitness_normalized)) %>%
  ungroup() %>%
  filter(ecology == "mutualism")  %>%
  lm(sd_fitness ~ E * M, data = .) %>% 
  summary()
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.02943    0.01070  96.242  < 2e-16 ***
#   ETRUE       -0.16279    0.01513 -10.762 9.80e-09 ***
#   MTRUE       -0.04007    0.01513  -2.649   0.0175 *  
#   ETRUE:MTRUE  0.14745    0.02139   6.893 3.62e-06 ***
fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  group_by(community, ecology, replicate, E, M, treatment) %>%
  summarize(mean_fitness = sd(fitness_normalized)) %>%
  ungroup() %>%
  filter(ecology == "competition") %>%
  lm(mean_fitness ~ E * M, data = .) %>% 
  summary()


#############
# Gene level linear regression analysis 
#############

# unhash to repeat analysis #

# analyze_fitness = function(fitness, response_var = "fitness_normalized"){
#   # does the two-factor lm on each gene and saves the data, then does BH corrections
#   # it does this on
#   results = data.frame()
#   fitness =  fitness  %>%
#     rename_(response_var = "response_var")
#   for (gene in unique(fitness$locusId)){
#     fitness %>%
#       filter(locusId == gene)%>%
#       lm(response_var ~ E + M + E:M, data = .) %>% # this is the equivalent of lm(response_var ~ E * M, data = .)
#       summary() %>%
#       coef %>%
#       data.frame() -> cf
#     names(cf) = c("Estimate", "SE", "t","p")
#     cf$term = row.names(cf)
#     cf$gene = gene
#     results = results %>%
#       rbind(cf)
#   }
#   # Benjamin-Hochberg corrections
#   results %>%
#     group_by(term) %>%
#     mutate(padj = p.adjust(p, method = "BH")) %>%
#     ungroup() %>%
#     filter(!is.na(padj))
# }
# # #
# results = analyze_fitness(fitness %>%
#                             filter(treatment %in% mutualism))
# write_csv(results, "data/mutualism_two_way_anova_stats.csv")

## competitive anova analysis

# results_comp = analyze_fitness(fitness %>%
#                             filter(treatment %in% competition))
#
# write_csv(results_comp, "data/competition_two_way_anova_stats.csv")


