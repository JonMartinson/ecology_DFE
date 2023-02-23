rm(list = ls())
library(tidyverse)
library(ggridges)
options(dplyr.summarise.inform = FALSE)

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
fitness %>%
  mutate(community = factor(community,
                            levels = rev(c("S","SE","SM","SEM")))) %>%
  ggplot(aes(x = fitness_normalized, y = community, fill = community))+
  geom_density_ridges2(rel_min_height = 0,
                       jittered_points = TRUE, panel_scaling = FALSE,
                       position = position_points_jitter(height = 0.05, yoffset = -0.1), 
                       point_alpha = 0.05, scale = 0.8, size = 0.3,
                       point_size = 0.2)+
  theme_bw(12)+
  scale_fill_discrete(guide = "none")+
  facet_wrap(~ecology)+
  labs(x = "gene-level fitness effect",
       y = "community members")
ggsave("plots/barseq_dfe_density_plus_raincloud.png",
       dpi = 300, width = 5, height = 2.5)

fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  group_by(community, ecology, locusId) %>%
  summarize(fitness_normalized = mean(fitness_normalized)) %>%
  ggplot(aes(x = fitness_normalized, fill = community))+
  geom_histogram(binwidth = 0.25)+
    theme_bw(12)+
  scale_y_continuous(limits = c(0, 100),
                     oob = scales::squish,
                     breaks = c(0,50,100),
                     labels = c("0", "50", "100+"))+
  scale_fill_discrete(guide = "none")+
  facet_grid(community~ecology)+
  labs(x = "gene-level fitness effect",
       y = "# genes")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("plots/barseq_dfe_squashed_histogram.png",
       dpi = 300, width = 5, height = 2.5)

fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  group_by(community, ecology, locusId) %>%
  summarize(fitness_normalized = mean(fitness_normalized)) %>%
  ggplot(aes(x = fitness_normalized, fill = community))+
  geom_histogram(binwidth = 0.25)+
  theme_bw(12)+
  scale_y_continuous(trans = "log2")+
  scale_fill_discrete(guide = "none")+
  facet_grid(community~ecology)+
  labs(x = "gene-level fitness effect",
       y = "# genes")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("plots/barseq_dfe_log2_histogram.png",
       dpi = 300, width = 5, height = 2.5)

fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  group_by(community, ecology, locusId) %>%
  summarize(fitness_normalized = mean(fitness_normalized)) %>%
  ggplot(aes(x = fitness_normalized, fill = community))+
  geom_histogram(binwidth = 0.25)+
  theme_bw(12)+
  scale_fill_discrete(guide = "none")+
  facet_grid(community~ecology)+
  labs(x = "gene-level fitness effect",
       y = "# genes")+
  scale_y_continuous(breaks = c(0, 1000, 2000))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave("./plots/barseq_dfe_histogram_unsquashed.png",
       dpi = 300, width = 5, height = 2.5)

fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  group_by(community, ecology, replicate) %>%
  summarize(mean_fitness = mean(fitness_normalized)) %>%
  ggplot(aes(x = community, y = mean_fitness, color = community))+
  stat_summary(fun = mean, shape = "-", size = 3, color = "black")+
  geom_point(shape = 1)+
  facet_wrap(~ecology)+
  theme_bw(12)+
  scale_color_discrete(guide = "none")+
  #geom_text(aes(label = replicate))+
  labs(x = "community members",
       y = "mean fitness effect")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(size = 0.5, color = "gray"))

ggsave("./plots/barseq_dfe_mean_fitness_effect_mean_stat.png",
       dpi = 300, width = 4, height = 2.5)

#statistics--parametric regression
fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  group_by(community, ecology, replicate, E, M, treatment) %>%
  summarize(mean_fitness = mean(fitness_normalized)) %>%
  ungroup() %>%
  filter(ecology == "mutualism") %>%
  lm(mean_fitness ~ E * M, data = .) %>%
  summary # 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.204542   0.003404 -60.083  < 2e-16 ***
#   ETRUE        0.041439   0.004814   8.607 2.12e-07 ***
#   MTRUE        0.029967   0.004814   6.224 1.22e-05 ***
#   ETRUE:MTRUE -0.049679   0.006809  -7.296 1.79e-06 ***
fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  group_by(community, ecology, replicate, E, M, treatment) %>%
  summarize(mean_fitness = mean(fitness_normalized)) %>%
  ungroup() %>%
  filter(ecology == "competition") %>%
  lm(mean_fitness ~ E + M, data = .) %>%
  summary
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.229037   0.003203 -71.503   <2e-16 ***
#   ETRUE       -0.005785   0.003769  -1.535    0.144    
# MTRUE       -0.001636   0.003769  -0.434    0.670  

#####Stdev in fitness effects
fitness %>%
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
       y = "stdev fitness effect")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(size = 0.5, color = "gray"))
ggsave("./plots/barseq_dfe_sd_fitness_effect_mean_stat.png",
       dpi = 300, width = 4, height = 2.5)

#statistics
fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  group_by(community, ecology, replicate, E, M, treatment) %>%
  summarize(mean_fitness = sd(fitness_normalized)) %>%
  ungroup() %>%
  filter(ecology == "mutualism")  %>%
  lm(mean_fitness ~ E * M, data = .) %>% 
  summary
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
  summary



#####coefficient of variation in fitness effects
fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  group_by(community, ecology, replicate) %>%
  summarize(mean_fitness = sd(fitness_normalized) / abs(mean(fitness_normalized))) %>%
  ggplot(aes(x = community, y = mean_fitness, color = community))+
  stat_summary(fun = mean, shape = "-", size = 3, color = "black")+
  geom_point(shape = 1)+
  facet_wrap(~ecology)+
  theme_bw(12)+
  scale_color_discrete(guide = "none")+
  #geom_text(aes(label = replicate))+
  labs(x = "community members",
       y = "CV of fitness effect")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(size = 0.5, color = "gray"))
ggsave("./plots/barseq_dfe_cv_fitness_effect_mean_stat.png",
       dpi = 300, width = 4, height = 2.5)

#statistics
fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  group_by(community, ecology, replicate, E, M, treatment) %>%
  summarize(mean_fitness = sd(fitness_normalized) / abs(mean(fitness_normalized))) %>%
  ungroup() %>%
  filter(ecology == "mutualism")  %>%
  lm(mean_fitness ~ E * M, data = .) %>% 
  summary
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.03384    0.06024  83.562  < 2e-16 ***
#   ETRUE        0.28299    0.08519   3.322  0.00432 ** 
#   MTRUE        0.64330    0.08519   7.551 1.16e-06 ***
#   ETRUE:MTRUE -0.63038    0.12048  -5.232 8.22e-05 ***

fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  group_by(community, ecology, replicate, E, M, treatment) %>%
  summarize(mean_fitness = sd(fitness_normalized) / mean(fitness_normalized)) %>%
  ungroup() %>%
  filter(ecology == "competition")   %>%
  lm(mean_fitness ~ E * M, data = .) %>% 
  summary
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -4.87686    0.05090 -95.813   <2e-16 ***
#   ETRUE        0.07933    0.07198   1.102    0.288    
# MTRUE       -0.09182    0.07198  -1.276    0.222    
# ETRUE:MTRUE  0.11096    0.10493   1.057    0.307  

#############
# Two-way ANOVA analysis
#############

#unhash to repeat analysis

# analyze_fitness = function(fitness, response_var = "fitness_normalized"){
#   # does the two-factor lm on each gene and saves the data, then does BH corrections
#   # it does this on
#   results = data.frame()
#   fitness =  fitness  %>%
#     rename_(response_var = "response_var")
#   for (gene in unique(fitness$locusId)){
#     fitness %>%
#       filter(locusId == gene)%>%
#       lm(response_var ~ E + M + E:M, data = .) %>%
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


##### mutualism anova analysis #####

results = read_csv("./data/S0_barseq_wetlab/mutualism_two_way_anova_stats.csv")

##############
# Why are the summary stats the way they are?
##############


# HOW MANY LOW-EFFECT GENES?? 
# try to make plots like I did for the t-test data, but using the two-way analysis
# to do this, I will be examining both significance as well as whether the final
# estimate is more extreme that -1 or 1. The main question right now is how to do it
# for the SEM treatment. At the moment, if no main effects besides intercept are significant,
# it just counts those signficant (low) intercepts. If any main effects, or the interact,
# are significant, then it examines the fitness estimate of the SEM community. An alternative
# would be to sum just the significant main effects (and intercept) if the interation
# is not significant. I don't think it would matter much either way. 
mutualism_low_fitness_genes = results %>%
  group_by(gene) %>%
  summarize(S_low = padj[term == "(Intercept)"] < 0.05 & 
              Estimate[term == "(Intercept)"] < -1,
            S_estimate = Estimate[term == "(Intercept)"],
            E_sig = padj[term == "ETRUE"] < 0.05,
            SE_estimate = sum(Estimate[term %in% c("(Intercept)", "ETRUE")]),
            M_sig = padj[term == "MTRUE"] < 0.05,
            SM_estimate = sum(Estimate[term %in% c("(Intercept)", "MTRUE")]),
            inter_sig = padj[term == "ETRUE:MTRUE"] < 0.05,
            SEM_estimate = sum(Estimate)) %>%
  ungroup() %>%
  summarize(S = sum(S_low),
            SE = sum(c(S_low[!E_sig], 
                       SE_estimate[S_low & E_sig] < -1,
                       SE_estimate[!S_low & E_sig] < -1)),
            SM = sum(c(S_low[!M_sig], 
                       SM_estimate[S_low & M_sig] < -1,
                       SM_estimate[!S_low & M_sig] < -1)),
            SEM = sum(c(S_low[!E_sig & !M_sig & !inter_sig],
                        SEM_estimate[E_sig | M_sig | inter_sig] < -1))) %>%
  pivot_longer(cols = -c(),
               names_to = "community", 
               values_to = "low_fitness_genes") %>%
  mutate(community = factor(community,
                            levels = c("S", "SE", "SM", "SEM"))) 
mutualism_low_fitness_genes %>%
  ggplot(aes(x = community, y = low_fitness_genes))+
  geom_bar(stat = "identity")+
  theme_bw(12)+
  labs(x = "community members",
       y = "low-fitness knockouts")
ggsave("./plots/barseq_twoway_low_fitness_knockouts_mutualism.png",
       dpi = 300, width = 2.5, height = 2.5)
mutualism_low_fitness_genes %>%
  select(low_fitness_genes) %>%
  pull %>%
  chisq.test() # p = 0.71

#rm(mutualism_low_fitness_genes)

# if we bin by low, moderate, and high fitness, is it always lower in monoculture?
fitness %>%
  group_by(community, locusId, rep, ecology) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  group_by(community, ecology, rep) %>%
  summarize(low = mean(fitness[fitness < -1]),
            neutral = mean(fitness[(fitness >= -1) & (fitness < 1)]),
            high = mean(fitness[fitness >= 1])) %>%
  pivot_longer(cols = c(low, neutral, high), names_to = "fitness_level",
               values_to = "bin_means") %>%
  mutate(fitness_level = ifelse(fitness_level == "low", "fitness < -1",
                                ifelse(fitness_level == "neutral",
                                       "|fitness| < 1", "fitness > 1"))) %>%
  mutate(fitness_level = factor(fitness_level,
                                levels = c("fitness < -1",
                                           "|fitness| < 1",
                                           "fitness > 1"))) %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  ggplot(aes(x = community, y = bin_means, color = community))+
  geom_point()+
  stat_summary(fun = mean, shape = "-", size = 2, color = "black")+
  facet_wrap(ecology~fitness_level, scales = "free_y")+
  theme_bw(12)+
  labs(y = "mean fitness in bin")
ggsave("./plots/barseq_dfe_mean_by_fitness_bin.png",
       dpi = 300, width = 6, height = 2.5)

# If we remove the low-fitness bin, do the mean fitness effect differences
# disappear?
fitness %>%
  filter(fitness_normalized > -1 ) %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  group_by(community, ecology, replicate) %>%
  summarize(mean_fitness = mean(fitness_normalized)) %>%
  ggplot(aes(x = community, y = mean_fitness, color = community))+
  stat_summary(fun = mean, shape = "-", size = 3, color = "black")+
  geom_point(shape = 1)+
  facet_wrap(~ecology)+
  theme_bw(12)+
  scale_color_discrete(guide = "none")+
  #geom_text(aes(label = replicate))+
  labs(x = "community members",
       y = "mean fitness effect")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(size = 0.5, color = "gray"))
ggsave("./plots/barseq_mean_fitness_effect_after_removing_low_fitness_kos.png",
       dpi = 300, width = 4, height = 2.5)
fitness %>%
  filter(fitness_normalized > -1 ) %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  group_by(community, ecology, replicate, E, M) %>%
  summarize(mean_fitness = mean(fitness_normalized)) %>%
  lm(mean_fitness ~ E + M, data = .) %>% 
  summary()

# another go at the same idea, using difference plots
fitness %>%
  group_by(community, locusId, ecology) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  filter(ecology == "mutualism") %>%
  pivot_wider(names_from = "community", values_from = "fitness") %>%
  ggplot(aes(x = S, y = SE))+
  geom_point(alpha = 0.5)+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw(12)+
  xlim(-8.5,2)+
  ylim(-8.5,2)
ggsave("./plots/barseq_scatter_S_vs_SE.png",
       dpi = 300, width = 3, height = 2.5)
fitness %>%
  group_by(community, locusId, ecology) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  filter(ecology == "mutualism") %>%
  pivot_wider(names_from = "community", values_from = "fitness") %>%
  ggplot(aes(x = S, y = SM))+
  geom_point(alpha = 0.5)+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw(12)+
  xlim(-8.5,2)+
  ylim(-8.5,2)
ggsave("./plots/barseq_scatter_S_vs_SM.png",
       dpi = 300, width = 3, height = 2.5)

fitness %>%
  group_by(community, locusId, ecology) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  filter(ecology == "mutualism") %>%
  pivot_wider(names_from = "community", values_from = "fitness") %>%
  summarize(lower = c(sum(SE < S), sum(S < SE)),
            lower_species = c("SE","S")) %>%
  ggplot(aes(x = lower_species, y = lower))+
  geom_bar(stat = "identity")+
  theme_bw(12)+
  labs(x = "community members",
       y = "# of knockouts\nwith lower fitness estimate")
ggsave("./plots/barseq_number_knockouts_lower_in_S_vs_SE_mutualism.png",
       dpi = 300, width = 2.5, height = 2.5)
fitness %>%
  group_by(community, locusId, ecology) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  filter(ecology == "mutualism") %>%
  pivot_wider(names_from = "community", values_from = "fitness") %>%
  summarize(lower = c(sum(SM < S), sum(S < SM)),
            lower_species = c("SM","S")) %>%
  ggplot(aes(x = lower_species, y = lower))+
  geom_bar(stat = "identity")+
  theme_bw(12)+
  labs(x = "community members",
       y = "# of knockouts\nwith lower fitness estimate")
ggsave("./plots/barseq_number_knockouts_lower_in_S_vs_SM_mutualism.png",
       dpi = 300, width = 2.5, height = 2.5)

dist_pt_to_line = function(x, y, a = -1, b = 1, c = 0){
  # calculates distance from a point (x,y) to a line
  # defined by ax + by + c = 0. Default is a y=x line.
  abs(a* x  + b * y + c) / sqrt(a ^2 + b ^2)
}
fitness %>%
  group_by(community, locusId, ecology) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  filter(ecology == "mutualism") %>%
  pivot_wider(names_from = "community", values_from = "fitness") %>%
  mutate(lower = ifelse(S < SE, "S", "SE")) %>%
  mutate(dist = dist_pt_to_line(S, SE)) %>%
  group_by(lower) %>%
  summarize(mean_dist = mean(dist),
            sem_dist = sd(dist) / sqrt(n()-1)) %>%
  ggplot(aes(x = lower, y = mean_dist,
             ymax = mean_dist + sem_dist,
             ymin = mean_dist - sem_dist))+
  geom_bar(stat = "identity", width = 0.67)+
  geom_errorbar(width = 0.5)+
  labs(x = "community where\nKO had lower fitness",
       y = "distance from y = x")+
  theme_bw(12)
ggsave("./plots/barseq_S_vs_SE_dist_from_y_x.png",
       dpi = 300, width = 2.5, height = 2.5)

fitness %>%
  group_by(community, locusId, ecology) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  filter(ecology == "mutualism") %>%
  pivot_wider(names_from = "community", values_from = "fitness") %>%
  mutate(lower = ifelse(S < SM, "S", "SM")) %>%
  mutate(dist = dist_pt_to_line(S, SM)) %>%
  group_by(lower) %>%
  summarize(mean_dist = mean(dist),
            sem_dist = sd(dist) / sqrt(n()-1)) %>%
  ggplot(aes(x = lower, y = mean_dist,
             ymax = mean_dist + sem_dist,
             ymin = mean_dist - sem_dist))+
  geom_bar(stat = "identity", width = 0.67)+
  geom_errorbar(width = 0.5)+
  labs(x = "community where\nKO had lower fitness",
       y = "distance from y = x")+
  theme_bw(12)
ggsave("./plots/barseq_S_vs_SM_dist_from_y_x.png",
       dpi = 300, width = 2.5, height = 2.5)

#### 
# The DFE parameters could change due to just one gene, or each community
# could have the same # of lowest genes, but differ in the magnitude.
# Here, I ask, how frequently is S monoculture the lowest fitness gene?
####
mutualism_lowest_genes = fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  filter(treatment %in% mutualism) %>%
  group_by(locusId, community) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  group_by(locusId) %>%
  summarize(lowest = community[fitness == min(fitness)]) 
mutualism_lowest_genes %>%
  ggplot(aes(x = lowest))+
  geom_bar()+
  theme_bw(12)+
  labs(x = "community members", 
       y = "# genes where community\nhad lowest fitness")
ggsave("./plots/barseq_count_of_times_community_had_lowest_fitness_mutualism.png",
       dpi = 300, width = 2.5, height = 2.5)

mutualism_lowest_genes %>%
  group_by(lowest) %>%
  summarize(n = n()) %>%
  select(n) %>%
  pull() %>%
  chisq.test() # p = 5.97e-10
rm(mutualism_lowest_genes)
############
# so S monoculture doesn't have the lowest fitness for more genes. Is its fitness
# much lower when it is low fitness?
############
fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  filter(treatment %in% mutualism) %>%
  group_by(locusId, community) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  group_by(locusId) %>%
  summarize(community = community[fitness == min(fitness)],
            fitness_diff = min(fitness) - mean(fitness[fitness != min(fitness)])) %>%
  ggplot(aes(x = fitness_diff, color = community))+
  geom_density(bw = 0.5, alpha = 0.25, size = 1)+
  scale_y_continuous(oob = scales::squish,
                     limits = c(0, 0.125),
                     breaks = c(0, 0.125), labels = c("0", "0.125+"))+
  theme_bw(12)+
  labs(x = "fitness difference",
       y = "density")
ggsave("./plots/barseq_dfe_mutualism_cost_when_lowest_density.png",
       dpi = 300, width = 3, height = 2.5)


## is there a difference in the mean fitness difference when we examine
## the cases where each community had the lowest fitness?

fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  filter(treatment %in% mutualism) %>%
  group_by(locusId, community) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  group_by(locusId) %>%
  summarize(lowest = community[fitness == min(fitness)],
            fitness_diff = min(fitness) - mean(fitness[fitness != min(fitness)])) %>%
  group_by(lowest) %>%
  summarize(d = mean(fitness_diff),
            dmax = mean(fitness_diff) + sd(fitness_diff) / sqrt(n()-1),
            dmin = mean(fitness_diff) - sd(fitness_diff) / sqrt(n()-1)) %>%
  ggplot(aes(x = lowest, y = d, ymax = dmax, ymin = dmin))+
  geom_errorbar()+
  geom_point()+
  theme_bw(12)+
  labs(x = "community members",
       y = "fitness difference")
ggsave("./plots/barseq_dfe_fitness_difference_when_lowest_mean_sem.png",
       dpi = 300, width = 2.5, height = 2.5)


fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  filter(treatment %in% mutualism) %>%
  group_by(locusId, community) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  group_by(locusId) %>%
  summarize(lowest = community[fitness == min(fitness)],
            fitness_diff = min(fitness) - mean(fitness[fitness != min(fitness)])) %>%
  ungroup() %>%
  aov(fitness_diff ~ lowest, data = .) %>% 
  TukeyHSD

############# 
# Repeat the above, but for high fitness, because fewer high-fitness
# knockouts could also explain the lower mean fitness effect in S mono
##############

# 1. parse out the high-fitness knockouts from the two-way data
mutualism_high_fitness_genes = results %>%
  group_by(gene) %>%
  summarize(S_high = padj[term == "(Intercept)"] < 0.05 & 
              Estimate[term == "(Intercept)"] > 1,
            S_estimate = Estimate[term == "(Intercept)"],
            E_sig = padj[term == "ETRUE"] < 0.05,
            SE_estimate = sum(Estimate[term %in% c("(Intercept)", "ETRUE")]),
            M_sig = padj[term == "MTRUE"] < 0.05,
            SM_estimate = sum(Estimate[term %in% c("(Intercept)", "MTRUE")]),
            inter_sig = padj[term == "ETRUE:MTRUE"] < 0.05,
            SEM_estimate = sum(Estimate)) %>%
  ungroup() %>%
  summarize(S = sum(S_high),
            SE = sum(c(S_high[!E_sig], 
                       SE_estimate[S_high & E_sig] > 1,
                       SE_estimate[!S_high & E_sig] > 1)),
            SM = sum(c(S_high[!M_sig], 
                       SM_estimate[S_high & M_sig] > 1,
                       SM_estimate[!S_high & M_sig] > 1)),
            SEM = sum(c(S_high[!E_sig & !M_sig & !inter_sig],
                        SEM_estimate[E_sig | M_sig | inter_sig] > 1))) %>%
  pivot_longer(cols = -c(),
               names_to = "community", 
               values_to = "high_fitness_genes") %>%
  mutate(community = factor(community,
                            levels = c("S", "SE", "SM", "SEM"))) 
mutualism_high_fitness_genes %>%
  ggplot(aes(x = community, y = high_fitness_genes))+
  geom_bar(stat = "identity")+
  theme_bw(12)+
  labs(x = "community members",
       y = "high-fitness knockouts")
ggsave("./plots/barseq_twoway_high_fitness_knockouts_mutualism.png",
       dpi = 300, width = 2.5, height = 2.5)
mutualism_high_fitness_genes %>%
  select(high_fitness_genes) %>%
  pull %>%
  chisq.test() # p = 0.98
#rm(mutualism_high_fitness_genes)


mutualism_highest_fitness = fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  filter(treatment %in% mutualism) %>%
  group_by(locusId, community) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  group_by(locusId) %>%
  summarize(highest = community[fitness == max(fitness)]) 
mutualism_highest_fitness %>%
  ggplot(aes(x = highest))+
  geom_bar()+
  theme_bw(12)+
  labs(x = "community members", 
       y = "# genes where community\nhad highest fitness")
ggsave("./plots/barseq_count_of_times_community_had_highest_fitness_mutualism.png",
       dpi = 300, width = 2.5, height = 2.5)
mutualism_highest_fitness %>%
  group_by(highest) %>%
  summarize(n = n()) %>%
  select(n) %>% pull() %>%
  chisq.test() # p = 0.0006
  
fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  filter(treatment %in% mutualism) %>%
  group_by(locusId, community) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  group_by(locusId) %>%
  summarize(community = community[fitness == max(fitness)],
            fitness_diff = max(fitness) -mean(fitness[fitness != max(fitness)]))%>%
  ggplot(aes(x = fitness_diff, color = community))+
  geom_density(bw = 0.5, alpha = 0.25, size = 1)+
  scale_y_continuous(oob = scales::squish,
                     limits = c(0, 0.125),
                     breaks = c(0, 0.125), labels = c("0", "0.125+"))+
  theme_bw(12)+
  labs(x = "fitness difference",
       y = "density")
ggsave("./plots/barseq_dfe_mutualism_benefit_when_highest_density.png",
       dpi = 300, width = 3, height = 2.5)


## is there a difference in the mean fitness difference when we examine
## the cases where each community had the lowest fitness?

fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  filter(treatment %in% mutualism) %>%
  group_by(locusId, community) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  group_by(locusId) %>%
  summarize(community = community[fitness == max(fitness)],
            fitness_diff = max(fitness) -mean(fitness[fitness != max(fitness)])) %>%
  group_by(community) %>%
  summarize(d = mean(fitness_diff),
            dmax = mean(fitness_diff) + sd(fitness_diff) / sqrt(n()-1),
            dmin = mean(fitness_diff) - sd(fitness_diff) / sqrt(n()-1)) %>%
  ggplot(aes(x = community, y = d, ymax = dmax, ymin = dmin))+
  geom_errorbar()+
  geom_point()+
  theme_bw(12)+
  labs(x = "community members",
       y = "fitness difference")
ggsave("./plots/barseq_dfe_fitness_difference_when_highest_mean_sem.png",
       dpi = 300, width = 2.5, height = 2.5)


fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  filter(treatment %in% mutualism) %>%
  group_by(locusId, community) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  group_by(locusId) %>%
  summarize(community = community[fitness == max(fitness)],
            fitness_diff = max(fitness) - mean(fitness[fitness != max(fitness)])) %>%
  ungroup() %>%
  aov(fitness_diff ~ community, data = .) %>% 
  TukeyHSD



fitness %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM"))) %>%
  filter(treatment %in% mutualism) %>%
  group_by(locusId, community) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  group_by(locusId) %>%
  summarize(lowest = community[fitness == max(fitness)],
            fitness_diff = max(fitness) -mean(fitness[fitness != max(fitness)]))%>%
  ggplot(aes(x = fitness_diff))+
  geom_histogram(binwidth = 0.0625)+
  scale_y_continuous(limits = c(0, 100),
                     oob = scales::squish,
                     breaks = c(0,50,100),
                     labels = c("0", "50", "100+"))+
  theme_bw(12)+
  facet_wrap(~lowest, ncol = 1)+
  labs(x = "fitness (when highest) - mean fitness")

###############
# A step further, what are the differences caused by the different 
# dependencies? Still at the distributional level, but past the point of summary stats
# More like, which genes are important in M, E?
#
# the order of figs here is 1) hierarchical clustering, 2) venn diagram of low fitness genes
# 3) genes that ameliorate, 4) fitness of genes that ameliorate, if any, genes that make things worse,
# 5) how many additive genes there are (2 sig main effects, no interaction), how many non-additive there are, colored by sign of interaction
###############

# hierarchical clustering
library(ggdendro)
mut_clust = fitness %>%
  filter(treatment %in% mutualism) %>%
  mutate(id = paste(community, rep)) %>%
  select(locusId, id, fitness_normalized) %>%
  pivot_wider(names_from = locusId, values_from = fitness_normalized) %>%
  column_to_rownames(var = "id") %>%
  dist(method = "euclidean") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram(hang = 0.1) %>%
  dendro_data()
ggplot()+
  geom_segment(data = segment(mut_clust),
               aes(x = x, y = y,
                   xend = xend, yend = yend))+
  geom_text(data = segment(mut_clust) %>%
              filter((x - floor(x)) == 0) %>%
              select(-y) %>%
              left_join(mut_clust$labels) %>%
              select(-y) %>% mutate(y = yend),
            aes(x = x, y = y, label = label), 
            angle = 90, hjust = 1)+
  ylim(10,60)+
  theme_void(12)
ggsave("./plots/barseq_mutualism_hier.png",
       dpi = 300, width = 4, height = 2.5)

# venn diagram
results %>%
  group_by(gene) %>%
  summarize(S_low = padj[term == "(Intercept)"] < 0.05 & 
              Estimate[term == "(Intercept)"] < -1,
            S_estimate = Estimate[term == "(Intercept)"],
            E_sig = padj[term == "ETRUE"] < 0.05,
            SE_estimate = sum(Estimate[term %in% c("(Intercept)", "ETRUE")]),
            M_sig = padj[term == "MTRUE"] < 0.05,
            SM_estimate = sum(Estimate[term %in% c("(Intercept)", "MTRUE")]),
            inter_sig = padj[term == "ETRUE:MTRUE"] < 0.05,
            SEM_estimate = sum(Estimate)) %>%
  group_by(gene) %>%
  summarize(S = S_low,
            SE = ifelse(!E_sig, S_low, SE_estimate < -1),
            SM = ifelse(!M_sig, S_low, SM_estimate < -1),
            SEM = ifelse(!E_sig & !M_sig & !inter_sig,
                         S_low, SEM_estimate < -1)) %>%
  data.frame() -> low_fitness_genes
low_list = list()
for (community in names(low_fitness_genes)[2:5]){
  low_list[[community]] = low_fitness_genes$gene[low_fitness_genes[,community]]
}
ggVennDiagram::ggVennDiagram(low_list,
                             label = "count",
                             edge_size = 1,
                             label_alpha = 0)+
  scale_fill_distiller(palette = "RdBu", direction = 1)+
  theme_bw(12)+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank())
ggsave("./plots/barseq_mutualism_lowfitness_venn_twoway.png",
       dpi = 300, width = 4, height = 2.5)

### number of low-fitness monoculture knockouts ameliorated
results %>%
  mutate(locusId = gene) %>%
  group_by(locusId) %>%
  filter(padj[term == "(Intercept)"] < 0.05 & 
           Estimate[term == "(Intercept)"] < 0)  %>%# now have cases where mono is < 0
  filter(any(padj[!(term %in% c("(Intercept)", "ETRUE:MTRUE"))] < 0.05))  %>% # and where at least one other term is sign
  filter(any(Estimate[!(term %in% c("(Intercept)", "ETRUE:MTRUE")) & padj < 0.05] > 0)) %>% # and where that effect is positve
  ungroup() %>%
  filter(term %in% c("ETRUE", "MTRUE")) %>%
  group_by(locusId) %>%
  summarize(effects = ifelse(all(padj < 0.05) & all(Estimate > 0), "both",
                             ifelse(padj[term == "ETRUE"] < 0.05 & Estimate[term == "ETRUE"] > 0,
                                    "only E","only M"))) %>% #write_csv("./data/low_monoculture_fitness_genes_rescued_by_mutualists.csv")
  ggplot(aes(x = effects))+
  geom_bar()+
  theme_bw(12)+
  labs(x = "species that ameliorate\nlow-fitness knockout",
       y = "genes ameliorated")
ggsave("./plots/barseq_mutualism_number_of_genes_ameliorated_by_E_M_or_both.png",
       dpi = 300, width = 2.5, height = 2.5)

### amount of amelioration
results %>%
  mutate(locusId = gene) %>%
  group_by(locusId) %>%
  filter(padj[term == "(Intercept)"] < 0.05 & 
           Estimate[term == "(Intercept)"] < 0)  %>%# now have cases where mono is < 0
  filter(any(padj[!(term %in% c("(Intercept)", "ETRUE:MTRUE"))] < 0.05))  %>% # and where at least one other term is sign
  filter(any(Estimate[!(term %in% c("(Intercept)", "ETRUE:MTRUE")) & padj < 0.05] > 0)) %>% # and where that effect is positve
  ungroup() %>%
  filter(term %in% c("ETRUE", "MTRUE")) %>%
  filter(padj < 0.05 & Estimate > 0) %>%
  group_by(locusId) %>%
  mutate(effects = ifelse(n() == 2, "both",
                          ifelse(term == "ETRUE", "only E", "only M"))) %>%
  mutate(benefactor = ifelse(term == "ETRUE", "E", "M")) %>%
  ggplot(aes(x = effects, y = Estimate, color = benefactor))+
  geom_point(size = 0.5,
             position = position_jitterdodge(dodge.width = 0.5,
                                             jitter.width = 0.1))+
  geom_boxplot(position = position_dodge(width = 0.5), width = 0.2,
               alpha = 0.7)+
  scale_color_brewer(type = "qual",
                     guide = guide_legend(""))+
  labs(x = "species that ameliorate\nlow-fitness knockout",
       y = "fitness improvement due\nto species' presence")+
  theme_bw(12)+
  theme(legend.position = c(0.75,0.85),
        legend.background = element_blank())
ggsave("./plots/barseq_mutualism_degree_of_amelioration_by_E_M_or_both.png",
       dpi = 300, width = 2.5, height = 2.5)


# genes with which are worsened by E, M
results %>%
  mutate(locusId = gene) %>%
  filter(term %in% c("ETRUE", "MTRUE")) %>%
  filter(padj < 0.05 & Estimate < 0) %>%
  group_by(locusId) %>%
  mutate(effects = ifelse(n() == 2, "both",
                             ifelse(term == "ETRUE", "only E","only M"))) %>%
  ggplot(aes(x = effects))+
  geom_bar()+
  theme_bw(12)+
  labs(x = "species that lower\nknockout fitness",
       y = "genes worsened")
ggsave("./plots/barseq_mutualism_number_of_genes_worsened_by_E_M_or_both.png",
       dpi = 300, width = 2.5, height = 2.5)


# amount of worsening
results %>%
  mutate(locusId = gene) %>%
  filter(term %in% c("ETRUE", "MTRUE")) %>%
  filter(padj < 0.05 & Estimate < 0) %>%
  group_by(locusId) %>%
  mutate(effects = ifelse(n() == 2, "both",
                          ifelse(term == "ETRUE", "only E","only M"))) %>%
  mutate(benefactor = ifelse(term == "ETRUE", "E", "M")) %>%
  ggplot(aes(x = effects, y = Estimate, color = benefactor))+
  geom_point(size = 0.5,
             position = position_jitterdodge(dodge.width = 0.5,
                                             jitter.width = 0.1))+
  geom_boxplot(position = position_dodge(width = 0.5), width = 0.2,
               alpha = 0.7)+
  scale_color_brewer(type = "qual",
                     guide = guide_legend(""))+
  labs(x = "species that lower\nknockout fitness",
       y = "fitness detriment due\nto species' presence")+
  theme_bw(12)+
  theme(legend.position = c(0.25,0.25),
        legend.background = element_blank())
ggsave("./plots/barseq_mutualism_degree_of_worsening_by_E_M_or_both.png",
       dpi = 300, width = 2.5, height = 2.5)

results %>%
  mutate(locusId = gene) %>%
  filter(term %in% c("ETRUE", "MTRUE")) %>%
  filter(padj < 0.05 & Estimate < 0) %>%
  group_by(locusId) %>%
  mutate(effects = ifelse(n() == 2, "both",
                          ifelse(term == "ETRUE", "only E","only M"))) %>%
  mutate(benefactor = ifelse(term == "ETRUE", "E", "M")) %>%
  t.test(Estimate ~ term, data = .) # p = 0.01994

## additivity
results %>%
  group_by(gene) %>%
  filter(padj[term == "ETRUE"] < 0.05 & padj[term == "MTRUE"] < 0.05) %>%
  ungroup() %>%
  filter(term == "ETRUE:MTRUE") %>%
  mutate(interaction = ifelse(padj >= 0.05,
                              "additive","non-additive")) %>%
  mutate(effect = ifelse(interaction == "additive", NA,
                         ifelse(Estimate > 0, "positive", "negative"))) %>%
  ggplot(aes(x = interaction, fill = effect))+
  scale_fill_discrete(breaks = c("negative", "positive"),
                      guide = guide_legend("change from\nadditivity"))+
  geom_bar()+
  labs(x = "interaction between E and M",
       y = "# of genes")+
  theme_bw(12)
  
##### MORE ADDITIVITY EXAMINATION
#### Aim for something like Jon's upset plot, but only 
#### for the genes where there are significant interactions
#### Actually, just make a flowchart for significance and effects.
#### show some small histograms for sign / magnitude of effects

# how many neutral?
results %>%
  group_by(gene) %>%
  filter(all(padj > 0.05)) %>%
  ungroup() %>% filter(term == "(Intercept)") %>% nrow # 3180
# leaving...
results %>%
  group_by(gene) %>%
  filter(!all(padj > 0.05)) %>%
  ungroup() %>% filter(term == "(Intercept)") %>% nrow # 370
# how many universal fitness effects? 
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] <= 0.05) &
           all(padj[term != "(Intercept)"] > 0.05)) %>%
  ungroup() %>% filter(term == "(Intercept)") %>% nrow # 235
# how many NOT intercept but with E, M, or E:M?
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] >= 0.05) &
           any(padj[term != "(Intercept)"] < 0.05)) %>%
  ungroup() %>% filter(term == "(Intercept)") %>% nrow # 18
# Not intercept, yes E
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] >= 0.05) &
           (padj[term == "ETRUE"] < 0.05) &
           (padj[term == "MTRUE"] >= 0.05 )) %>%
  ungroup() %>% filter(term == "ETRUE")
# Not intercept, yes E, yes interaction
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] >= 0.05) &
           (padj[term == "ETRUE"] < 0.05) &
           (padj[term == "MTRUE"] >= 0.05 ) &
           (padj[term == "ETRUE:MTRUE"] < 0.05)) %>%
  ungroup() %>% filter(term == "ETRUE")
# Not intercept, yes M
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] >= 0.05) &
           (padj[term == "MTRUE"] < 0.05) &
           (padj[term == "ETRUE"] >= 0.05 )) %>%
  ungroup() %>% filter(term == "MTRUE")
# Not intercept, yes M, yes interaction
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] >= 0.05) &
           (padj[term == "MTRUE"] < 0.05) &
           (padj[term == "ETRUE"] >= 0.05 ) &
           (padj[term == "ETRUE:MTRUE"] < 0.05)) %>%
  ungroup() %>% filter(term == "ETRUE:MTRUE")
# Not intercept, yes E and M
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] >= 0.05) &
           (padj[term == "ETRUE:MTRUE"] < 0.05)  ) %>%
  ungroup() 
# yes intercept, yes E (not M)?
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] < 0.05) &
           (padj[term == "ETRUE"] < 0.05) &
           (padj[term == "MTRUE"] >= 0.05 )) %>%
  ungroup() %>% filter(term == "ETRUE")  %>% nrow # 23
# yes intercept, yes E (not M), yes interaction??
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] < 0.05) &
           (padj[term == "ETRUE"] < 0.05) &
           (padj[term == "MTRUE"] >= 0.05 ) &
           (padj[term == "ETRUE:MTRUE"] < 0.05)) %>%
  ungroup() %>% filter(term == "ETRUE:MTRUE")  
# yes intercept, yes M (not E)?
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] < 0.05) &
           (padj[term == "MTRUE"] < 0.05) &
           (padj[term == "ETRUE"] >= 0.05 )) %>%
  ungroup() %>% filter(term == "MTRUE")  %>% nrow # 23
#yes intercept, yes M, no E, yes interaction
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] < 0.05) &
           (padj[term == "ETRUE"] >= 0.05) &
           (padj[term == "MTRUE"] < 0.05 ) &
           (padj[term == "ETRUE:MTRUE"] < 0.05)) %>%
  ungroup() %>% filter(term == "ETRUE:MTRUE")  
#yes intercept, yes M AND E
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] < 0.05) &
           (padj[term == "ETRUE"] < 0.05) &
           (padj[term == "MTRUE"] < 0.05 )) %>%
  ungroup() %>% filter(term == "ETRUE:MTRUE") 
#yes intercept, yes M AND E AND interaction
results %>%
  group_by(gene) %>%
  filter(all(padj < 0.05)) %>%
  ungroup() %>% filter(term == "ETRUE:MTRUE") %>%
  summarize(reduced = sum(Estimate < 0),
            increased = sum(Estimate >= 0))
#yes intercept, no M no E yes interaction
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] < 0.05) &
           (padj[term == "ETRUE"] >= 0.05) &
           (padj[term == "MTRUE"] >= 0.05 ) &
           (padj[term == "ETRUE:MTRUE"] < 0.05))

###PLOTS--line plots
# 1. universal effects (sign. intercept only)
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] < 0.05) &
  all(padj[term != "(Intercept)"] >= 0.05)) %>%
  ungroup() %>%
  filter(term == "(Intercept)") %>%
  ggplot(aes(x = Estimate, y = 0))+
  geom_point(shape = "|", size =3)+
  theme_minimal(12)+
  scale_y_continuous(breaks = c())+
  geom_hline(yintercept = 0)+
  scale_x_continuous(breaks = c(-8, -4, 0))+
  labs(x = "", y = "")
ggsave("./plots/barseq_hist_mutualism_universal_fitnesses.png",
       dpi = 300, width = 1.5, height = 0.75)
# Not intercept, yes E
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] >= 0.05) &
           (padj[term == "ETRUE"] < 0.05) &
           (padj[term == "MTRUE"] >= 0.05 )) %>%
  ungroup() %>% filter(term == "ETRUE")%>%
  ggplot(aes(x = Estimate, y = 0))+
  geom_point(shape = "|", size =3)+
  theme_minimal(12)+
  scale_y_continuous(breaks = c())+
  geom_hline(yintercept = 0)+
  scale_x_continuous(limits = c(-2, 0), breaks = c(-2, -1, 0))+
  labs(x = "", y = "")
ggsave("./plots/barseq_hist_mutualism_E_only_fitnesses.png",
       dpi = 300, width = 1.5, height = 0.75)
# Not intercept, yes M
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] >= 0.05) &
           (padj[term == "MTRUE"] < 0.05) &
           (padj[term == "ETRUE"] >= 0.05 )) %>%
  ungroup() %>% filter(term == "MTRUE")%>%
  ggplot(aes(x = Estimate, y = 0))+
  geom_point(shape = "|", size =3)+
  theme_minimal(12)+
  scale_y_continuous(breaks = c())+
  geom_hline(yintercept = 0)+
  scale_x_continuous(limits = c(-4, 2), breaks = c(-4, -2, 0))+
  labs(x = "", y = "")
ggsave("./plots/barseq_hist_mutualism_M_only_fitnesses.png",
       dpi = 300, width = 1.5, height = 0.75)
# Yes intercept, yes E or M or RM
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] < 0.05) &
          any(padj[term != "(Intercept)"] < 0.05)) %>%
  ungroup() %>% filter(term == "(Intercept)")%>%
  ggplot(aes(x = Estimate, y = 0))+
  geom_point(shape = "|", size =3)+
  theme_minimal(12)+
  scale_y_continuous(breaks = c())+
  geom_hline(yintercept = 0)+
  scale_x_continuous(limits = c(-8, 2), breaks = c(-8, -4, 0))+
  labs(x = "", y = "")
ggsave("./plots/barseq_hist_mutualism_intercept_when_atleastoneothersig.png",
       dpi = 300, width = 1.5, height = 0.75)
# yes intercept, yes E, no M
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] < 0.05) &
          (padj[term == "ETRUE"] < 0.05) &
           (padj[term == "MTRUE"] >= 0.05)) %>%
  ungroup() %>% filter(term %in% c("ETRUE", "(Intercept)"))%>%
  select(-c(SE,t,p,padj)) %>%
  pivot_wider(names_from = term, values_from = Estimate) %>%
  ggplot(aes(x = ETRUE, y = 0, color = `(Intercept)` < 0))+
  geom_point(shape = "|", size =3)+
  scale_color_manual(values = c("purple", "green"),
                     guide = "none")+
  theme_minimal(12)+
  scale_y_continuous(breaks = c())+
  geom_hline(yintercept = 0)+
  scale_x_continuous(limits = c(-2, 3), breaks = c(-2, 0, 2))+
  labs(x = "", y = "")
ggsave("./plots/barseq_hist_mutualism_interceptyes_Eyes_Mno_Eestimate.png",
       dpi = 300, width = 1.5, height = 0.75)
# yes intercept, yes M, no E
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] < 0.05) &
           (padj[term == "MTRUE"] < 0.05) &
           (padj[term == "ETRUE"] >= 0.05)) %>%
  ungroup() %>% filter(term %in% c("MTRUE", "(Intercept)"))%>%
  select(-c(SE,t,p,padj)) %>%
  pivot_wider(names_from = term, values_from = Estimate) %>%
  ggplot(aes(x = MTRUE, y = 0, color = `(Intercept)` < 0))+
  geom_point(shape = "|", size =3)+
  scale_color_manual(values = c("purple", "green"),
                     guide = "none")+
  theme_minimal(12)+
  scale_y_continuous(breaks = c())+
  geom_hline(yintercept = 0)+
  scale_x_continuous(limits = c(-6, 3), breaks = c(-2, 0, 2))+
  labs(x = "", y = "")
ggsave("./plots/barseq_hist_mutualism_interceptyes_Myes_Eno_Mestimate.png",
       dpi = 300, width = 1.5, height = 0.75)
# yes intercept, yes M, yes E
results %>%
  group_by(gene) %>%
  filter((padj[term == "(Intercept)"] < 0.05) &
           (padj[term == "MTRUE"] < 0.05) &
           (padj[term == "ETRUE"] < 0.05)) %>%
  ungroup() %>% filter(term %in% c("MTRUE","ETRUE", "(Intercept)"))%>%
  select(-c(SE,t,p,padj)) %>%
  pivot_wider(names_from = term, values_from = Estimate) %>%
  pivot_longer(c(MTRUE,ETRUE), names_to = "term",
               values_to = "Estimate") %>%
  mutate(y = ifelse(term == "MTRUE", 0, 1)) %>%
  ggplot(aes(x = Estimate, y = y, color = `(Intercept)` < 0))+
  geom_point(shape = "|", size =3)+
  scale_color_manual(values = c("purple", "green"),
                     guide = "none")+
  theme_minimal(12)+
  scale_y_continuous(breaks = c())+
  geom_hline(yintercept = c(0,1))+
  scale_x_continuous(limits = c(-3,5), breaks = c(-2, 0, 2))+
  labs(x = "", y = "")
ggsave("./plots/barseq_hist_mutualism_interceptyes_Myes_Eyes_Mestimate.png",
       dpi = 300, width = 1.5, height = 1)

