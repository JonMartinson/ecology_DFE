
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
  annotate("text", x = c(-6, -1), y = c(1, -7.5), label = c("n = 57 y > x", "n = 28 y < x"))+
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
  labs(x = "fitness in S",y = "fitness in SE")+
  theme(legend.position = "none",
        panel.grid = element_blank())
ggsave("./plots/barseq_S_vs_SE_scatter_with_sig_color.png",
       dpi = 300, width = 3, height = 2.5)


fitness %>%
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
  annotate("text", x = c(-6, -1), y = c(1, -7.5), label = c("n = 79 y > x", "n = 17 y < x"))+
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
  labs(x = "fitness in S",y = "fitness in SM")+
  theme(legend.position = "none",
        panel.grid = element_blank())
ggsave("./plots/barseq_S_vs_SM_scatter_with_sig_color.png",
       dpi = 300, width = 3, height = 2.5)

fitness %>%
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
  mutate(designation = ifelse(padj_SE < 0.05 & padj_SM < 0.05, "both",
                              ifelse(padj_SE < 0.05, "E",
                                     ifelse(padj_SM < 0.05, "M", "none")))) %>% 
  mutate(designation = factor(designation,
                              levels = c("none","E","M", "both"))) %>%
  ggplot(aes(x = fitness_mean_S,
             xmin = fitness_mean_S - fitness_sd_S,
             xmax = fitness_mean_S + fitness_sd_S,
             y = fitness_mean_SEM,
             ymax = fitness_mean_SEM + fitness_sd_SEM,
             ymin = fitness_mean_SEM - fitness_sd_SEM,
             color = designation, shape = padj_SEM < 0.05))+
  scale_color_manual(values = rev(c("gray", "green", "blue", "red")),
                     guide = guide_legend(title = "significant effect of: "))+
  scale_shape_discrete(guide = "none")+
  geom_errorbar(width = 0, color = "gray")+
  geom_errorbarh(height = 0, color = "gray")+
  geom_point(data = .%>% filter(padj_SE > 0.05 & padj_SM > 0.05 & padj_SEM > 0.05))+
  geom_point(data = .%>% filter(padj_SE <= 0.05 | padj_SM <= 0.05 | padj_SEM <= 0.05))+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw(12)+
  scale_x_continuous(expand = c(0,0), limits = c(-8, 1.5),
                     oob = scales::squish)+
  scale_y_continuous(expand = c(0,0), limits = c(-8, 1.5),
                     oob = scales::squish)+
  labs(x = "fitness in S",y = "fitness in SEM")+
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")
ggsave("./plots/barseq_S_vs_SEM_scatter_with_sig_color.png",
       dpi = 300, width = 3, height = 3)

fitness %>%
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
  mutate(designation = ifelse(padj_SE < 0.05 & padj_SM < 0.05, "both",
                              ifelse(padj_SE < 0.05, "E",
                                     ifelse(padj_SM < 0.05, "M", "none")))) %>% 
  mutate(designation = factor(designation,
                              levels = c("none","E","M", "both"))) %>%
  ggplot(aes(x = fitness_mean_SE,
             xmin = fitness_mean_SE - fitness_sd_SE,
             xmax = fitness_mean_SE + fitness_sd_SE,
             y = fitness_mean_SEM,
             ymax = fitness_mean_SEM + fitness_sd_SEM,
             ymin = fitness_mean_SEM - fitness_sd_SEM,
             color = designation))+
  scale_color_manual(values = c("none" = "darkgray", 
                                "E" = "green", 
                                "M" = "blue", "both" = "red"),
                     guide = guide_legend(title = "effect of: "))+  
  stat_poly_eq(mapping = use_label("R2"), color = 1)+
  geom_errorbar(width = 0, color = "lightgray")+
  geom_errorbarh(height = 0, color = "lightgray")+
  geom_point(data = .%>%filter(designation == "none"))+
  geom_point(data = .%>%filter(designation != "none"))+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw(12)+
  scale_x_continuous(expand = c(0,0), limits = c(-8, 1.5),
                     oob = scales::squish)+
  scale_y_continuous(expand = c(0,0), limits = c(-8, 1.5),
                     oob = scales::squish)+
  labs(x = "fitness in SE",y = "fitness in SEM")+
  theme(legend.position = "top",
        legend.spacing.x = unit(0, "mm"),
        panel.grid = element_blank())
ggsave("./plots/barseq_SE_vs_SEM_scatter_with_r2.png",
       dpi = 300, width = 3, height = 3)

fitness %>%
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
  mutate(designation = ifelse(padj_SE < 0.05 & padj_SM < 0.05, "both",
                              ifelse(padj_SE < 0.05, "E",
                                     ifelse(padj_SM < 0.05, "M", "none")))) %>% 
  mutate(designation = factor(designation,
                              levels = c("none","E","M", "both"))) %>%
  ggplot(aes(x = fitness_mean_SM,
             xmin = fitness_mean_SM - fitness_sd_SM,
             xmax = fitness_mean_SM + fitness_sd_SM,
             y = fitness_mean_SEM,
             ymax = fitness_mean_SEM + fitness_sd_SEM,
             ymin = fitness_mean_SEM - fitness_sd_SEM,
             color = designation))+
  scale_color_manual(values = c("none" = "darkgray", 
                                "E" = "green", 
                                "M" = "blue", "both" = "red"),
                     guide = guide_legend(title = "effect of: "))+  
  stat_poly_eq(mapping = use_label("R2"), color = 1)+
  geom_errorbar(width = 0, color = "lightgray")+
  geom_errorbarh(height = 0, color = "lightgray")+
  geom_point(data = .%>%filter(designation == "none"))+
  geom_point(data = .%>%filter(designation != "none"))+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw(12)+
  scale_x_continuous(expand = c(0,0), limits = c(-8, 1.5),
                     oob = scales::squish)+
  scale_y_continuous(expand = c(0,0), limits = c(-8, 1.5),
                     oob = scales::squish)+
  labs(x = "fitness in SM",y = "fitness in SEM")+
  theme(legend.position = "top",
        legend.spacing.x = unit(0, "mm"),
        panel.grid = element_blank())
ggsave("./plots/barseq_SM_vs_SEM_scatter_with_r2.png",
       dpi = 300, width = 3, height = 3)

fitness %>%
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
  mutate(designation = ifelse(padj_SE < 0.05 & padj_SM < 0.05, "both",
                              ifelse(padj_SE < 0.05, "E",
                                     ifelse(padj_SM < 0.05, "M", "none")))) %>% 
  mutate(designation = factor(designation,
                              levels = c("none","E","M", "both"))) %>%
  ggplot(aes(x = fitness_mean_SE,
             xmin = fitness_mean_SE - fitness_sd_SE,
             xmax = fitness_mean_SE + fitness_sd_SE,
             y = fitness_mean_SM,
             ymax = fitness_mean_SM + fitness_sd_SM,
             ymin = fitness_mean_SM - fitness_sd_SM,
             color = designation))+
  scale_color_manual(values = c("none" = "darkgray", 
                                "E" = "green", 
                                "M" = "blue", "both" = "red"),
                     guide = guide_legend(title = "effect of: "))+  
  stat_poly_eq(mapping = use_label("R2"), color = 1)+
  geom_errorbar(width = 0, color = "lightgray")+
  geom_errorbarh(height = 0, color = "lightgray")+
  geom_point(data = .%>%filter(designation == "none"))+
  geom_point(data = .%>%filter(designation != "none"))+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw(12)+
  scale_x_continuous(expand = c(0,0), limits = c(-8, 1.5),
                     oob = scales::squish)+
  scale_y_continuous(expand = c(0,0), limits = c(-8, 1.5),
                     oob = scales::squish)+
  labs(x = "fitness in SE",y = "fitness in SM")+
  theme(legend.position = "top",
        legend.spacing.x = unit(0, "mm"),
        panel.grid = element_blank())
ggsave("./plots/barseq_SE_vs_SM_scatter_with_r2.png",
       dpi = 300, width = 3, height = 3)


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
         SEM_predicted = Estimate_S + Estimate_SE + Estimate_SM) %>% # n_SM = 1758, n_S = 1792
  mutate(designation = ifelse(padj_SE < 0.05 & padj_SM < 0.05, "both",
                              ifelse(padj_SE < 0.05, "E",
                                     ifelse(padj_SM < 0.05, "M", "none")))) %>% 
  mutate(designation = factor(designation,
                              levels = c("none","E","M", "both"))) %>%
  filter(designation != "none") %>%
  ggplot(aes(x = SEM_predicted,
             y = fitness_mean_SEM,
             color = padj_SEM < 0.05))+ 
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "purple"),
                     guide = guide_legend("non-additive:"))+
  scale_shape_discrete(guide = "none")+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw(12)+
  scale_x_continuous(limits = c(-8, 1.5),
                     oob = scales::squish)+
  scale_y_continuous(limits = c(-8, 1.5),
                     oob = scales::squish)+
  labs(x = "SEM predicted",y = "SEM actual")+
  theme(legend.position = "top",
        legend.spacing.x = unit(0, "mm"),
        panel.grid = element_blank())
ggsave("./plots/barseq_SEM_predicted_vs_actual_only_sign_effects.png",
       dpi = 300, width = 3, height = 3)


fitness %>%
  filter(gene_name == "mdh") %>%
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
  ggplot(aes(x = community, y = fitness_mean,
             ymax = fitness_mean + fitness_sd,
             ymin = fitness_mean - fitness_sd))+
  scale_x_discrete(limits = c("S","SE","SM","SEM"))+
  geom_point()+
  geom_errorbar(width = 0.3)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_bw(12)+
  labs(y = "relative fitness",
       title = "mdh gene, mutualism")



### competition analysis
# change stats to competitive anova
stats = read_csv("./data/competition_two_way_anova_stats.csv")

stats = stats %>% 
    mutate(locusId = gene) %>%
    left_join(gene_names)

# competitive S vs SE


fitness %>%
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
    annotate("text", x = c(-6, -1), y = c(1, -7.5), label = c("n = 25 y > x", "n = 50 y < x"))+ #these values are correct for comp
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
    labs(x = "fitness in S",y = "fitness in SE")+
    theme(legend.position = "none",
          panel.grid = element_blank())

ggsave("./plots/COMPETITION_barseq_S_vs_SE_scatter_with_sig_color.png",
       dpi = 300, width = 3, height = 2.5)


#S vs SM

# THIS CODE IS MODIFIED TO MATCH COLORS -- DON'T RERUN ON OTHER DATA
fitness %>%
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
           worsen = sum(padj_SM < 0.05 & (fitness_mean_SM < fitness_mean_S)))  %>% # ameliorate -- 79; worsen -- 17
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
    annotate("text", x = c(-6, -1), y = c(1, -7.5), label = c("n = 0 y > x", "n = 0 y < x"))+ # nothing significantly different
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
    labs(x = "fitness in S",y = "fitness in SM")+
    theme(legend.position = "none",
          panel.grid = element_blank())
ggsave("./plots/COMPETITION_barseq_S_vs_SM_scatter_with_sig_color.png",
       dpi = 300, width = 3, height = 2.5)

fitness %>%
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
    mutate(designation = ifelse(padj_SE < 0.05 & padj_SM < 0.05, "both",
                                ifelse(padj_SE < 0.05, "E",
                                       ifelse(padj_SM < 0.05, "M", "none")))) %>% 
    mutate(designation = factor(designation,
                                levels = c("none","E","M", "both"))) %>%
    ggplot(aes(x = fitness_mean_S,
               xmin = fitness_mean_S - fitness_sd_S,
               xmax = fitness_mean_S + fitness_sd_S,
               y = fitness_mean_SEM,
               ymax = fitness_mean_SEM + fitness_sd_SEM,
               ymin = fitness_mean_SEM - fitness_sd_SEM,
               color = designation, shape = padj_SEM < 0.05))+
    scale_color_manual(values = rev(c("gray", "green", "blue", "red")),
                       guide = guide_legend(title = "significant effect of: "))+
    scale_shape_discrete(guide = "none")+
    geom_errorbar(width = 0, color = "gray")+
    geom_errorbarh(height = 0, color = "gray")+
    geom_point(data = .%>% filter(padj_SE > 0.05 & padj_SM > 0.05 & padj_SEM > 0.05))+
    geom_point(data = .%>% filter(padj_SE <= 0.05 | padj_SM <= 0.05 | padj_SEM <= 0.05))+
    geom_abline(slope = 1, intercept = 0)+
    theme_bw(12)+
    scale_x_continuous(expand = c(0,0), limits = c(-8, 1.5),
                       oob = scales::squish)+
    scale_y_continuous(expand = c(0,0), limits = c(-8, 1.5),
                       oob = scales::squish)+
    labs(x = "fitness in S",y = "fitness in SEM")+
    theme(panel.grid.minor = element_blank(),
          legend.position = "top")


