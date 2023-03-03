
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

A <- fitness %>%
  filter(ecology == "mutualism") %>%
  group_by(locusId, E, M, community, gene_name) %>%
  summarize(fitness_mean = mean(fitness_normalized ),
            fitness_sd = sd(fitness_normalized ) / sqrt(4)) %>%
  left_join(stats %>%
              select(term, padj, locusId, gene_name) %>%
              mutate(community = ifelse(term == "(Intercept)",
                                        "S",
                                        ifelse(term == "ETRUE",
                                               "SE",
                                               ifelse(term == "MTRUE", "SM", "SEM"))))) %>%
  ungroup() %>%
  select(community, locusId, gene_name, fitness_mean, fitness_sd, padj) %>%
  pivot_wider(names_from = "community", values_from = c("fitness_mean", "fitness_sd", "padj")) %>%
  mutate(amelioration = sum(padj_SE < 0.05 & (fitness_mean_SE > fitness_mean_S)),
            worsen = sum(padj_SE < 0.05 & (fitness_mean_SE < fitness_mean_S))) %>% # ameliorate -- 57; worsen -- 28
  ggplot(aes(x = fitness_mean_S,
             xmin = fitness_mean_S - fitness_sd_S,
             xmax = fitness_mean_S + fitness_sd_S,
             y = fitness_mean_SE,
             ymax = fitness_mean_SE + fitness_sd_SE,
             ymin = fitness_mean_SE - fitness_sd_SE,
             color = padj_SE < 0.05))+
  scale_color_brewer(type = "qual", palette = 6, direction = -1)+
  geom_errorbar(width = 0, color = "gray")+
  geom_errorbarh(height = 0, color = "gray")+
  annotate("text", x = c(-6, -1), y = c(1, -7.5), label = c("57", "28"))+
  geom_point(data = .%>%
               filter(padj_SE >= 0.05))+
  geom_point(data = .%>%
               filter(padj_SE < 0.05))+
  geom_abline(slope = 1, intercept = 0)+
  #stat_poly_eq(mapping = use_label("R2"), color = 1)+
  
  theme_bw(12)+
  scale_x_continuous(expand = c(0,0), limits = c(-8, 1.5),
                     oob = scales::squish)+
  scale_y_continuous(expand = c(0,0), limits = c(-8, 1.5),
                     oob = scales::squish)+
  labs(x = "Fitness in monoculture",y = "Fitness in coculture with E")+
  theme(legend.position = "none",
        panel.grid = element_blank())



B <- fitness %>%
  filter(ecology == "mutualism") %>%
  group_by(locusId, E, M, community, gene_name) %>%
  summarize(fitness_mean = mean(fitness_normalized ),
            fitness_sd = sd(fitness_normalized ) / sqrt(4)) %>%
  left_join(stats %>%
              select(term, padj, locusId, gene_name) %>%
              mutate(community = ifelse(term == "(Intercept)",
                                        "S",
                                        ifelse(term == "ETRUE",
                                               "SE",
                                               ifelse(term == "MTRUE", "SM", "SEM"))))) %>%
  ungroup() %>%
  select(community, locusId, gene_name, fitness_mean, fitness_sd, padj) %>%
  pivot_wider(names_from = "community", values_from = c("fitness_mean", "fitness_sd", "padj")) %>%
  mutate(amelioration = sum(padj_SM < 0.05 & (fitness_mean_SM > fitness_mean_S)),
         worsen = sum(padj_SM < 0.05 & (fitness_mean_SM < fitness_mean_S))) %>% # ameliorate -- 79; worsen -- 17
  ggplot(aes(x = fitness_mean_S,
             xmin = fitness_mean_S - fitness_sd_S,
             xmax = fitness_mean_S + fitness_sd_S,
             y = fitness_mean_SM,
             ymax = fitness_mean_SM + fitness_sd_SM,
             ymin = fitness_mean_SM - fitness_sd_SM,
             color = padj_SM < 0.05))+
  scale_color_brewer(type = "qual", palette = 6, direction = -1)+
  #stat_poly_eq(mapping = use_label("R2"), color = 1)+
  
  geom_errorbar(width = 0, color = "gray")+
  geom_errorbarh(height = 0, color = "gray")+
  annotate("text", x = c(-6, -1), y = c(1, -7.5), label = c("79", "17"))+
  geom_point(data = .%>%
               filter(padj_SM >= 0.05))+
  geom_point(data = .%>%
               filter(padj_SM < 0.05))+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw(12)+
  scale_x_continuous(expand = c(0,0), limits = c(-8, 1.5),
                     oob = scales::squish)+
  scale_y_continuous(expand = c(0,0), limits = c(-8, 1.5),
                     oob = scales::squish)+
  labs(x = "Fitness in monoculture",y = "Fitness in coculture with M")+
  theme(legend.position = "none",
        panel.grid = element_blank())


#### competition analysis ####
# change stats to competitive anova
stats = read_csv("./data/competition_two_way_anova_stats.csv")

stats = stats %>% 
    mutate(locusId = gene) %>%
    left_join(gene_names)

# competitive S vs SE


C <- fitness %>%
    filter(ecology == "competition") %>%
    group_by(locusId, E, M, community, gene_name) %>%
    summarize(fitness_mean = mean(fitness_normalized ),
              fitness_sd = sd(fitness_normalized ) / sqrt(4)) %>%
    left_join(stats %>%
                  select(term, padj, locusId, gene_name) %>%
                  mutate(community = ifelse(term == "(Intercept)",
                                            "S",
                                            ifelse(term == "ETRUE",
                                                   "SE",
                                                   ifelse(term == "MTRUE", "SM", "SEM"))))) %>%
    ungroup() %>%
    select(community, locusId, gene_name, fitness_mean, fitness_sd, padj) %>%
    pivot_wider(names_from = "community", values_from = c("fitness_mean", "fitness_sd", "padj")) %>%
    mutate(amelioration = sum(padj_SE < 0.05 & (fitness_mean_SE > fitness_mean_S)),
           worsen = sum(padj_SE < 0.05 & (fitness_mean_SE < fitness_mean_S))) %>% # ameliorate -- 25; worsen -- 50
    ggplot(aes(x = fitness_mean_S,
               xmin = fitness_mean_S - fitness_sd_S,
               xmax = fitness_mean_S + fitness_sd_S,
               y = fitness_mean_SE,
               ymax = fitness_mean_SE + fitness_sd_SE,
               ymin = fitness_mean_SE - fitness_sd_SE,
               color = padj_SE < 0.05))+
    scale_color_brewer(type = "qual", palette = 6, direction = -1)+
    geom_errorbar(width = 0, color = "gray")+
    geom_errorbarh(height = 0, color = "gray")+
    annotate("text", x = c(-6, -1), y = c(1, -7.5), label = c("25", "50"))+ #these values are correct for comp
    geom_point(data = .%>%
                   filter(padj_SE >= 0.05))+
    geom_point(data = .%>%
                   filter(padj_SE < 0.05))+
    geom_abline(slope = 1, intercept = 0)+
    #stat_poly_eq(mapping = use_label("R2"), color = 1)+
    
    theme_bw(12)+
    scale_x_continuous(expand = c(0,0), limits = c(-8, 1.5),
                       oob = scales::squish)+
    scale_y_continuous(expand = c(0,0), limits = c(-8, 1.5),
                       oob = scales::squish)+
    labs(x = "Fitness in monoculture",y = "Fitness in coculture with E")+
    theme(legend.position = "none",
          panel.grid = element_blank())



#S vs SM

# THIS CODE IS MODIFIED TO MATCH COLORS -- DON'T RERUN ON OTHER DATA
D <- fitness %>%
    filter(ecology == "competition") %>%
    group_by(locusId, E, M, community, gene_name) %>%
    summarize(fitness_mean = mean(fitness_normalized ),
              fitness_sd = sd(fitness_normalized ) / sqrt(4)) %>%
    left_join(stats %>%
                  select(term, padj, locusId, gene_name) %>%
                  mutate(community = ifelse(term == "(Intercept)",
                                            "S",
                                            ifelse(term == "ETRUE",
                                                   "SE",
                                                   ifelse(term == "MTRUE", "SM", "SEM"))))) %>%
    ungroup() %>%
    select(community, locusId, gene_name, fitness_mean, fitness_sd, padj) %>%
    pivot_wider(names_from = "community", values_from = c("fitness_mean", "fitness_sd", "padj")) %>%
    mutate(amelioration = sum(padj_SM < 0.05 & (fitness_mean_SM > fitness_mean_S)),
           worsen = sum(padj_SM < 0.05 & (fitness_mean_SM < fitness_mean_S)))  %>% # 
    ggplot(aes(x = fitness_mean_S,
               xmin = fitness_mean_S - fitness_sd_S,
               xmax = fitness_mean_S + fitness_sd_S,
               y = fitness_mean_SM,
               ymax = fitness_mean_SM + fitness_sd_SM,
               ymin = fitness_mean_SM - fitness_sd_SM,
               color = locusId == 'STM0005'))+ 
    scale_color_brewer(type = "qual", palette = 6, direction = -1)+
    #stat_poly_eq(mapping = use_label("R2"), color = 1)+
    
    geom_errorbar(width = 0, color = "gray")+
    geom_errorbarh(height = 0, color = "gray")+
    annotate("text", x = c(-6, -1), y = c(1, -7.5), label = c("0", "0"))+ # nothing significantly different
    geom_point(data = .%>%
                   filter(padj_SM >= 0.05))+
    geom_point(data = .%>%
                   filter(padj_SM < 0.05))+
    geom_abline(slope = 1, intercept = 0)+
    theme_bw(12)+
    scale_x_continuous(expand = c(0,0), limits = c(-8, 1.5),
                       oob = scales::squish)+
    scale_y_continuous(expand = c(0,0), limits = c(-8, 1.5),
                       oob = scales::squish)+
    labs(x = "Fitness in monoculture",y = "Fitness in coculture with M")+
    theme(legend.position = "none",
          panel.grid = element_blank())



ggpubr::ggarrange(A,B,C,D, labels = c('A','B','C','D'))
ggsave("./plots/all_interaction_scatter_with_sig_color.png",
       dpi = 300, width = 6, height = 6)
