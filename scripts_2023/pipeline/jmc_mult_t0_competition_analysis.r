rm(list = ls())
fitness = read_csv("./data/S0_barseq_wetlab/S0_barseq_fitnesses_20220722.csv")
gene_names = read_csv("./data/external_metadata/uniprot_s_enterica_gene_names.csv") %>%
  mutate(multi_name = `Gene Names`,
         gene_name = `Gene Names (primary)`) %>%
  mutate(locusId = sapply(multi_name,
                           FUN = function(x){
                             str_split(x, " ")[[1]][length(str_split(x, " ")[[1]])]
                           })) %>%
  select(gene_name, locusId)


comp = fitness %>%
  filter(str_detect(treatment, "comp")) %>%
  group_by(treatment, locusId, rep, replicate, desc) %>%
  summarize(fitness_normalized = first(fitness_normalized)) %>%
  mutate(E = ifelse(str_detect(treatment, "E"), TRUE, FALSE), # add variable whether E or M is present in the treatment
         M = ifelse(str_detect(treatment, "M "), TRUE, FALSE)) 



comp %>%
  group_by(treatment, replicate,locusId, E, M) %>%
  summarize(fitness_normalized = mean(fitness_normalized)) %>%
  ggplot(aes(x = E, y = fitness_normalized,  fill = M))+
  geom_violin()+
  theme_bw(12)+
  scale_fill_discrete(guide = guide_legend("M. extorquens\npresent"))+
  labs(x = "E. coli present",
       y = "gene KO fitness")

ggsave("./plots/dfe_violin_plot_competition.png",
       dpi = 300, width = 4, height = 3)


comp %>%
  group_by(treatment, replicate, E, M) %>%
  summarize(mean_fitness = mean(fitness_normalized)) %>%
  ggplot(aes(x = E, y = mean_fitness,  shape = M,
             color = M))+
  geom_point()+
  geom_line(data = comp %>%
              group_by(treatment, replicate, E, M) %>%
              summarize(mean_fitness = mean(fitness_normalized)) %>%
              group_by(E, M) %>%
              summarize(mean_fitness = mean(mean_fitness)),
            aes(group = M))+
  theme_bw(12)+
  scale_color_discrete(guide = "none")+
  scale_shape_discrete(guide = guide_legend("M. extorquens\npresent"))+
  labs(x = "E. coli present",
       y = "mean fitness effect")

ggsave("./plots/dfe_means_competition.png",
       dpi = 300, width = 4, height = 3)

comp %>%
  group_by(treatment, replicate, E, M) %>%
  summarize(mean_fitness = mean(fitness_normalized)) %>%
  aov(mean_fitness ~ E * M, data = .) %>%
  summary()


comp %>%
  group_by(treatment, replicate, E, M) %>%
  summarize(sd_fitness = sd(fitness_normalized)) %>%
  ggplot(aes(x = E, y = sd_fitness,  shape = M,
             color = M))+
  geom_point()+
  geom_line(data = comp %>%
              group_by(treatment, replicate, E, M) %>%
              summarize(sd_fitness = sd(fitness_normalized)) %>%
              group_by(E, M) %>%
              summarize(sd_fitness = mean(sd_fitness)),
            aes(group = M))+
  theme_bw(12)+
  scale_color_discrete(guide = "none")+
  scale_shape_discrete(guide = guide_legend("M. extorquens\npresent"))+
  labs(x = "E. coli present",
       y = "SD fitness effect")

ggsave("./plots/dfe_sds_competition.png",
       dpi = 300, width = 4, height = 3)
comp %>%
  group_by(treatment, replicate, E, M) %>%
  summarize(sd_fitness = sd(fitness_normalized)) %>%
  aov(sd_fitness ~ E * M, data = .) %>%
  summary()


analyze_fitness_by_trt = function(fitness){
  results = data.frame()
  for (gene in unique(fitness$locusId)){
    for (trt in unique(fitness$treatment)){
      this_test = fitness %>%
        ungroup() %>%
        filter(locusId == gene,
               treatment == trt) %>%
        select(fitness_normalized) %>%
        unlist() %>%
        t.test() 

      this_result = data.frame(gene, treatment = trt,
                               estimate = unname(this_test$estimate),
                               conf_bottom = this_test$conf.int[1],
                               conf_top = this_test$conf.int[1],
                               p = this_test$p.value)
      results = rbind(results, this_result)
    }
  }
  results$p.adj = p.adjust(results$p, method = "BH")
  results
}
comp_results = analyze_fitness_by_trt(comp)

comp_results %>%
  mutate(E = ifelse(str_detect(treatment, "E"), TRUE, FALSE), # add variable whether E or M is present in the treatment
         M = ifelse(str_detect(treatment, "M "), TRUE, FALSE)) %>%
  filter(estimate < -1 & p.adj < 0.05) %>%
  ggplot(aes(x = E, fill = M))+
  geom_bar(position = position_dodge())+
  labs(x = "E. coli presence",
       y = "# low-fitness KOs")+
  scale_fill_discrete(guide = guide_legend(title = "presence of\nM. extorquenz"))+
  theme_bw(12)
ggsave("./plots/dfe_lowfitnessgenes_competition.png",
       dpi = 300, width = 4, height = 3)

comp_results %>%
  mutate(E = ifelse(str_detect(treatment, "E"), TRUE, FALSE), # add variable whether E or M is present in the treatment
         M = ifelse(str_detect(treatment, "M "), TRUE, FALSE)) %>%
  filter(estimate < -1 & p.adj < 0.05) %>%
  group_by(E, M) %>%
  summarize(low_fitness = n()) %>%
  ungroup() %>%
  select(low_fitness) %>%
  unlist() %>%
  chisq.test(simulate.p.value = TRUE, B = 1e5) # X-squared = 3, df = NA, p-value = 0.3944

comp_results %>%
  mutate(E = ifelse(str_detect(treatment, "E"), TRUE, FALSE), # add variable whether E or M is present in the treatment
         M = ifelse(str_detect(treatment, "M "), TRUE, FALSE)) %>%
  filter(estimate > 1 & p.adj < 0.05) %>%
  ggplot(aes(x = E, fill = M))+
  geom_bar(position = position_dodge())+
  labs(x = "E. coli presence",
       y = "# high-fitness KOs")+
  scale_fill_discrete(guide = guide_legend(title = "presence of\nM. extorquenz"))+
  theme_bw(12)
ggsave("./plots/dfe_highfitnessgenes_competition.png",
       dpi = 300, width = 4, height = 3)


comp_results %>%
  mutate(E = ifelse(str_detect(treatment, "E"), TRUE, FALSE), # add variable whether E or M is present in the treatment
         M = ifelse(str_detect(treatment, "M "), TRUE, FALSE)) %>%
  filter(estimate > 1 & p.adj < 0.05) %>%
  group_by(E, M) %>%
  summarize(low_fitness = n()) %>%
  ungroup() %>%
  select(low_fitness) %>%
  unlist() %>%
  chisq.test(simulate.p.value = TRUE, B = 1e5) # X-squared = 8.5556, df = NA, p-value = 0.03546
