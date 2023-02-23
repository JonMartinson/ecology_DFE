
rm(list = ls())
library(tidyverse)
library(ggpmisc)

load("./data/fitness.Rdat")
load("./data/gene_names.Rdat")
stats = read_csv("./data/mutualism_two_way_anova_stats.csv")
gene_names$locusId[gene_names$gene_name == "ilvA"] = "STM3905"
gene_names$locusId[gene_names$gene_name == "ilvC"] = "STM3909"
gene_names$locusId[gene_names$gene_name == "ilvD"] = "STM3904"
gene_names$locusId[gene_names$gene_name == "ilvE"] = "STM3903"

gene_names = gene_names |>
  rbind(data.frame(gene_name = c("fadA","fadB"),
                   locusId = c("STM3982", "STM3983")))

fitness = fitness %>% left_join(gene_names)
stats = stats %>% 
  mutate(locusId = gene) %>%
  left_join(gene_names)


fitness %>%
  filter(ecology == "mutualism") %>%
  group_by(locusId, E, M, community, gene_name) %>%
  summarize(fitness_mean = mean(fitness_normalized ),
            fitness_sd = sd(fitness_normalized ) / sqrt(4)) %>%
  left_join(stats %>%
              select(term, padj, locusId, gene_name, Estimate) %>%
              mutate(community = ifelse(term == "(Intercept)",
                                        "S",
                                        ifelse(term == "ETRUE",
                                               "SE",
                                               ifelse(term == "MTRUE", "SM", "SEM"))))) %>%
  ungroup() %>%
  select(community, locusId, gene_name, fitness_mean, fitness_sd, Estimate, padj) %>%
  pivot_wider(names_from = "community", values_from = c("fitness_mean", "fitness_sd", "padj", "Estimate")) %>%
  mutate(n_SEM = sum(fitness_mean_SEM > fitness_mean_S),
         n_S = sum(fitness_mean_SEM < fitness_mean_S),
         additive = Estimate_S + Estimate_SE + Estimate_SM,
         average = Estimate_S + (Estimate_SE + Estimate_SM)/2,
         strongest = Estimate_S + ifelse(abs(Estimate_SE) > abs(Estimate_SM),
                                         Estimate_SE, Estimate_SM)) %>% # n_SM = 1758, n_S = 1792
  pivot_longer(cols = c(additive, average, strongest), names_to = "model",
               values_to = "prediction") %>%
  ggplot(aes(x = prediction,
             y = fitness_mean_SEM))+ 
  geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw(12)+
  stat_poly_eq(mapping = use_label("R2"), color = 1)+
  stat_smooth(method = "lm")+
    facet_wrap(~factor(model,levels = c("additive", "strongest", "average")))+
  labs(x = "SEM predicted",y = "SEM actual")+
  theme(legend.position = "top",
        legend.spacing.x = unit(0, "mm"),
        panel.grid = element_blank())
ggsave("./plots/mutualism_barseq_scatters_SEM_pred_vs_actual_diff_models.png",
       dpi = 300, height = 2.5, width = 8)

fitness %>%
  filter(ecology == "mutualism") %>%
  group_by(locusId, E, M, community, gene_name) %>%
  summarize(fitness_mean = mean(fitness_normalized ),
            fitness_sd = sd(fitness_normalized ) / sqrt(4)) %>%
  left_join(stats %>%
              select(term, padj, locusId, gene_name, Estimate) %>%
              mutate(community = ifelse(term == "(Intercept)",
                                        "S",
                                        ifelse(term == "ETRUE",
                                               "SE",
                                               ifelse(term == "MTRUE", "SM", "SEM"))))) %>%
  ungroup() %>%
  select(community, locusId, gene_name, fitness_mean, fitness_sd, Estimate, padj) %>%
  pivot_wider(names_from = "community", values_from = c("fitness_mean", "fitness_sd", "padj", "Estimate")) %>%
  mutate(additive = Estimate_S + Estimate_SE + Estimate_SM,
         average = Estimate_S + (Estimate_SE + Estimate_SM)/2,
         strongest = Estimate_S + ifelse(abs(Estimate_SE) > abs(Estimate_SM),
                                         Estimate_SE, Estimate_SM)) %>% # n_SM = 1758, n_S = 1792
  pivot_longer(cols = c(additive, average, strongest), names_to = "model",
               values_to = "prediction") %>%
  ggplot(aes(x = prediction,
             y = fitness_mean_SEM - prediction))+ 
  geom_point()+
  theme_bw(12)+
  geom_hline(yintercept = 0)+
  stat_smooth(method = "lm")+
  facet_wrap(~model)+
  labs(x = "SEM predicted",y = "SEM actual - predicted")+
  theme(legend.position = "top",
        legend.spacing.x = unit(0, "mm"),
        panel.grid = element_blank())+
  geom_text(data = data.frame(x = c(-5,-5,-5), y = c(-4,-4,-4),
                              model = c("additive", "average", "strongest"),
                              label = paste("RMSE =",c("0.401", "0.275", "0.323"))),
            aes(x=x,y=y,label=label))
ggsave("./plots/barseq_scatters_SEM_pred_diff_vs_actual_diff_models.png",
       dpi = 300, height = 2.5, width = 8)

fitness %>%
  filter(ecology == "mutualism") %>%
  group_by(locusId, E, M, community, gene_name) %>%
  summarize(fitness_mean = mean(fitness_normalized ),
            fitness_sd = sd(fitness_normalized ) / sqrt(4)) %>%
  left_join(stats %>%
              select(term, padj, locusId, gene_name, Estimate) %>%
              mutate(community = ifelse(term == "(Intercept)",
                                        "S",
                                        ifelse(term == "ETRUE",
                                               "SE",
                                               ifelse(term == "MTRUE", "SM", "SEM"))))) %>%
  ungroup() %>%
  select(community, locusId, gene_name, fitness_mean, fitness_sd, Estimate, padj) %>%
  pivot_wider(names_from = "community", values_from = c("fitness_mean", "fitness_sd", "padj", "Estimate")) %>%
  mutate(additive = Estimate_S + Estimate_SE + Estimate_SM,
         average = Estimate_S + (Estimate_SE + Estimate_SM)/2,
         strongest = Estimate_S + ifelse(abs(Estimate_SE) > abs(Estimate_SM),
                                         Estimate_SE, Estimate_SM)) %>% # n_SM = 1758, n_S = 1792
  pivot_longer(cols = c(additive, average, strongest), names_to = "model",
               values_to = "prediction") %>%
  group_by(model) %>%
  summarize(RMSE = sqrt(mean((fitness_mean_SEM - prediction)^2)))




####
# competitive #

stats = read_csv("./data/competition_two_way_anova_stats.csv")

stats = stats %>% 
    mutate(locusId = gene) %>%
    left_join(gene_names)

fitness %>%
    filter(ecology == "competition") %>%
    group_by(locusId, E, M, community, gene_name) %>%
    summarize(fitness_mean = mean(fitness_normalized ),
              fitness_sd = sd(fitness_normalized ) / sqrt(4)) %>%
    left_join(stats %>%
                  select(term, padj, locusId, gene_name, Estimate) %>%
                  mutate(community = ifelse(term == "(Intercept)",
                                            "S",
                                            ifelse(term == "ETRUE",
                                                   "SE",
                                                   ifelse(term == "MTRUE", "SM", "SEM"))))) %>%
    ungroup() %>%
    select(community, locusId, gene_name, fitness_mean, fitness_sd, Estimate, padj) %>%
    pivot_wider(names_from = "community", values_from = c("fitness_mean", "fitness_sd", "padj", "Estimate")) %>%
    mutate(n_SEM = sum(fitness_mean_SEM > fitness_mean_S),
           n_S = sum(fitness_mean_SEM < fitness_mean_S),
           additive = Estimate_S + Estimate_SE + Estimate_SM,
           average = Estimate_S + (Estimate_SE + Estimate_SM)/2,
           strongest = Estimate_S + ifelse(abs(Estimate_SE) > abs(Estimate_SM),
                                           Estimate_SE, Estimate_SM)) %>% # n_SM = 1758, n_S = 1792
    pivot_longer(cols = c(additive, average, strongest), names_to = "model",
                 values_to = "prediction") %>%
    ggplot(aes(x = prediction,
               y = fitness_mean_SEM))+ 
    geom_point()+
    geom_abline(slope = 1, intercept = 0)+
    theme_bw(12)+
    stat_poly_eq(mapping = use_label("R2"), color = 1)+
    stat_smooth(method = "lm")+
    facet_wrap(~factor(model,levels = c("additive", "strongest", "average")))+
    labs(x = "SEM predicted",y = "SEM actual")+
    theme(legend.position = "top",
          legend.spacing.x = unit(0, "mm"),
          panel.grid = element_blank())
ggsave("./plots/competition_barseq_scatters_SEM_pred_vs_actual_diff_models.png",
       dpi = 300, height = 2.5, width = 8)

