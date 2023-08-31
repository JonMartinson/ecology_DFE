rm(list = ls())
fitness = read_csv("data/S0_barseq_fitnesses_20220722.csv")
gene_names = read_csv("data/uniprot_s_enterica_gene_names.csv") %>%
  mutate(multi_name = `Gene Names`,
         gene_name = `Gene Names (primary)`) %>%
  mutate(locusId = sapply(multi_name,
                           FUN = function(x){
                             str_split(x, " ")[[1]][length(str_split(x, " ")[[1]])]
                           })) %>%
  select(gene_name, locusId)

early_mut <- c("T2 S mut", "T3 SE mut", "T5 SM mut", "T5 SEM mut") # early mutualism samples 

fitness = fitness %>%
  filter(treatment %in% early_mut) %>%
  group_by(treatment, locusId, rep, replicate, desc) %>%
  summarize(fitness_normalized = first(fitness_normalized)) %>%
  mutate(E = ifelse(str_detect(treatment, "E"), TRUE, FALSE), # add variable whether E or M is present in the treatment
         M = ifelse(str_detect(treatment, "M "), TRUE, FALSE)) 



fitness %>%
  group_by(treatment, replicate,locusId, E, M) %>%
  summarize(fitness_normalized = mean(fitness_normalized)) %>%
  ggplot(aes(x = E, y = fitness_normalized,  fill = M))+
  geom_violin()+
  theme_bw(12)+
  scale_fill_discrete(guide = guide_legend("M. extorquens\npresent"))+
  labs(x = "E. coli present",
       y = "gene KO fitness")

ggsave("./plots/dfe_violin_plot_mutualism.png",
       dpi = 300, width = 4, height = 3)


fitness %>%
  group_by(treatment, replicate, E, M) %>%
  summarize(mean_fitness = mean(fitness_normalized)) %>%
  ggplot(aes(x = E, y = mean_fitness,  shape = M,
             color = M))+
  geom_point()+
  geom_line(data = fitness %>%
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

ggsave("./plots/dfe_means_mutualism.png",
       dpi = 300, width = 4, height = 3)
fitness %>%
  group_by(treatment, replicate, E, M) %>%
  summarize(mean_fitness = mean(fitness_normalized)) %>%
  aov(mean_fitness ~ E * M, data = .) %>%
  summary()


fitness %>%
  group_by(treatment, replicate, E, M) %>%
  summarize(sd_fitness = sd(fitness_normalized)) %>%
  ggplot(aes(x = E, y = sd_fitness,  shape = M,
             color = M))+
  geom_point()+
  geom_line(data = fitness %>%
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

ggsave("./plots/dfe_sds_mutualism.png",
       dpi = 300, width = 4, height = 3)
fitness %>%
  group_by(treatment, replicate, E, M) %>%
  summarize(sd_fitness = sd(fitness_normalized)) %>%
  aov(sd_fitness ~ E * M, data = .) %>%
  summary()

# this one does t.tests
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
mut_results = analyze_fitness_by_trt(fitness)


mut_results %>%
  mutate(E = ifelse(str_detect(treatment, "E"), TRUE, FALSE), # add variable whether E or M is present in the treatment
         M = ifelse(str_detect(treatment, "M "), TRUE, FALSE)) %>%
  filter(estimate < -1 & p.adj < 0.05) %>%
  ggplot(aes(x = E, fill = M))+
  geom_bar(position = position_dodge())+
  labs(x = "E. coli presence",
       y = "# low-fitness KOs")+
  scale_fill_discrete(guide = guide_legend(title = "presence of\nM. extorquenz"))+
  theme_bw(12)
ggsave("./plots/dfe_lowfitnessgenes_mutualism.png",
       dpi = 300, width = 4, height = 3)

mut_results %>%
  mutate(E = ifelse(str_detect(treatment, "E"), TRUE, FALSE), # add variable whether E or M is present in the treatment
         M = ifelse(str_detect(treatment, "M "), TRUE, FALSE)) %>%
  filter(estimate < -1 & p.adj < 0.05) %>%
  group_by(E, M) %>%
  summarize(low_fitness = n()) %>%
  ungroup() %>%
  select(low_fitness) %>%
  unlist() %>%
  chisq.test(simulate.p.value = TRUE, B = 1e5) # X-squared = 7.7038, df = NA, p-value = 0.05205

mut_results %>%
  mutate(E = ifelse(str_detect(treatment, "E"), TRUE, FALSE), # add variable whether E or M is present in the treatment
         M = ifelse(str_detect(treatment, "M "), TRUE, FALSE)) %>%
  filter(estimate > 1 & p.adj < 0.05) %>%
  ggplot(aes(x = E, fill = M))+
  geom_bar(position = position_dodge())+
  labs(x = "E. coli presence",
       y = "# high-fitness KOs")+
  scale_fill_discrete(guide = guide_legend(title = "presence of\nM. extorquenz"))+
  theme_bw(12)
ggsave("./plots/dfe_highfitnessgenes_mutualism.png",
       dpi = 300, width = 4, height = 3)

mut_results %>%
  mutate(E = ifelse(str_detect(treatment, "E"), TRUE, FALSE), # add variable whether E or M is present in the treatment
         M = ifelse(str_detect(treatment, "M "), TRUE, FALSE)) %>%
  filter(estimate > 1 & p.adj < 0.05) %>%
  group_by(E, M) %>%
  summarize(low_fitness = n()) %>%
  ungroup() %>%
  select(low_fitness) %>%
  unlist() %>%
  chisq.test(simulate.p.value = TRUE, B = 1e5) # X-squared = 27.538, df = NA, p-value = 6e-05



analyze_fitness = function(fitness, response_var = "fitness_normalized"){
  # does the two-factor lm on each gene and saves the data, then does BH corrections
  # it does this on 
  results = data.frame()
  fitness =  fitness  %>%
    rename_(response_var = "response_var") 
  for (gene in unique(fitness$locusId)){
    fitness %>%
      filter(locusId == gene)%>%
      lm(response_var ~ E + M + E:M, data = .) %>%
      summary() %>%
      coef %>%
      data.frame() -> cf
    names(cf) = c("Estimate", "SE", "t","p")
    cf$term = row.names(cf)
    cf$gene = gene
    results = results %>%
      rbind(cf)
  }
  # Benjamin-Hochberg corrections
  results %>%
    group_by(term) %>%
    mutate(padj = p.adjust(p, method = "BH")) %>%
    ungroup() %>%
    filter(!is.na(padj))
}

results = analyze_fitness(fitness)



(a = fitness %>%
  mutate(treatment = ifelse(!E & !M, "monoculture",
                            ifelse(!E & M, "SM",
                                   ifelse(E & !M, "SE", "SEM")))) %>%
   group_by(locusId, treatment) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  pivot_wider(names_from = "treatment", values_from = "fitness") %>%
  left_join(results %>%
              mutate(locusId = gene) %>%
              select(term, locusId, padj) %>%
              pivot_wider(names_from = "term", values_from = "padj") %>%
              rename(intercept = `(Intercept)`)) %>%
  mutate(sign = ifelse(intercept < 0.05 & ETRUE >= 0.05, "intercept",
                ifelse(intercept < 0.05 & ETRUE < 0.05, "intercept and main E effect",
                ifelse(intercept >= 0.05 & ETRUE < 0.05, "main E effect", "other")))) %>%
  left_join(gene_names) %>%
  ggplot(aes(x = monoculture, y = SE, color = sign))+
  scale_color_discrete(guide =  guide_legend(title = ""))+
  geom_abline(intercept = 0, slope = 1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_point()+
  ggrepel::geom_text_repel(size = 2.5, seed = 1, force = 0.5, force_pull = 1.5,
                           color = "black", segment.color = "gray",
                           aes(label = gene_name))+
  theme_bw(12)+
  theme(legend.position = "top"))
ggsave("./plots/dfe_scatter_mono_v_SE_mutualism.png", plot = a,
       dpi = 300, width = 6, height = 6)


(b = fitness %>%
    mutate(treatment = ifelse(!E & !M, "monoculture",
                              ifelse(!E & M, "SM",
                                     ifelse(E & !M, "SE", "SEM")))) %>%
    group_by(locusId, treatment) %>%
    summarize(fitness = mean(fitness_normalized)) %>%
    pivot_wider(names_from = "treatment", values_from = "fitness") %>%
    left_join(results %>%
                mutate(locusId = gene) %>%
                select(term, locusId, padj) %>%
                pivot_wider(names_from = "term", values_from = "padj") %>%
                rename(intercept = `(Intercept)`)) %>%
    mutate(sign = ifelse(intercept < 0.05 & MTRUE >= 0.05, "intercept",
                         ifelse(intercept < 0.05 & MTRUE < 0.05, "intercept and main M effect",
                                ifelse(intercept >= 0.05 & MTRUE < 0.05, "main M effect", "other")))) %>%
    left_join(gene_names) %>%
    ggplot(aes(x = monoculture, y = SM, color = sign))+
    scale_color_discrete(guide =  guide_legend(title = ""))+
    geom_abline(intercept = 0, slope = 1)+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_point()+
    ggrepel::geom_text_repel(size = 2.5, seed = 1, force = 0.5, force_pull = 1.5,
                             color = "black", segment.color = "gray",
                             aes(label = gene_name))+
    theme_bw(12)+
    theme(legend.position = "top"))
ggsave("./plots/dfe_scatter_mono_v_SM_mutualism.png", plot = b,
       dpi = 300, width = 6, height = 6)

(c = fitness %>%
    mutate(treatment = ifelse(!E & !M, "monoculture",
                              ifelse(!E & M, "SM",
                                     ifelse(E & !M, "SE", "SEM")))) %>%
    group_by(locusId, treatment) %>%
    summarize(fitness = mean(fitness_normalized)) %>%
    pivot_wider(names_from = "treatment", values_from = "fitness") %>%
    left_join(results %>%
                mutate(locusId = gene) %>%
                select(term, locusId, padj) %>%
                pivot_wider(names_from = "term", values_from = "padj") %>%
                rename(intercept = `(Intercept)`)) %>%
    mutate(sign = ifelse(intercept < 0.05 & ETRUE >= 0.05 & MTRUE >= 0.05, "int.",
                  ifelse(intercept < 0.05 & ETRUE >= 0.05 & MTRUE < 0.05, "int. and main M",
                  ifelse(intercept < 0.05 & ETRUE < 0.05 & MTRUE >= 0.05, "int. and main E",
                  ifelse(intercept < 0.05 & ETRUE < 0.05 & MTRUE < 0.05, "int. and both", "other"))))) %>%
    left_join(gene_names) %>%
    ggplot(aes(x = monoculture, y = SEM, color = sign))+
    scale_color_discrete(guide =  guide_legend(title = ""))+
    geom_abline(intercept = 0, slope = 1)+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_point()+
    ggrepel::geom_text_repel(size = 2.5, seed = 1, force = 0.5, force_pull = 1.5,
                             color = "black", segment.color = "gray",
                             aes(label = gene_name))+
    theme_bw(12)+
    theme(legend.position = "top"))
ggsave("./plots/dfe_scatter_mono_v_SEM_mutualism.png", plot = c,
       dpi = 300, width = 6, height = 6)

(d = fitness %>%
    mutate(treatment = ifelse(!E & !M, "monoculture",
                              ifelse(!E & M, "SM",
                                     ifelse(E & !M, "SE", "SEM")))) %>%
    group_by(locusId, treatment) %>%
    summarize(fitness = mean(fitness_normalized)) %>%
    pivot_wider(names_from = "treatment", values_from = "fitness") %>%
    left_join(results %>%
                mutate(locusId = gene) %>%
                select(term, locusId, padj) %>%
                pivot_wider(names_from = "term", values_from = "padj") %>%
                rename(intercept = `(Intercept)`)) %>%
    mutate(sign = ifelse(intercept < 0.05 & ETRUE >= 0.05 & MTRUE >= 0.05, "int.",
                 ifelse(intercept < 0.05 & ETRUE >= 0.05 & MTRUE < 0.05, "int. and main M",
                  ifelse(intercept < 0.05 & ETRUE < 0.05 & MTRUE >= 0.05, "int. and main E",
                 ifelse(intercept < 0.05 & ETRUE < 0.05 & MTRUE < 0.05, "int. and both", "other"))))) %>%
    left_join(gene_names) %>%
    ggplot(aes(x = SE, y = SM, color = sign))+
    scale_color_discrete(guide =  guide_legend(title = ""))+
    geom_abline(intercept = 0, slope = 1)+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_point()+
    ggrepel::geom_text_repel(size = 2.5, seed = 1, force = 0.5, force_pull = 1.5,
                             color = "black", segment.color = "gray",
                             aes(label = gene_name))+
    theme_bw(12)+
    theme(legend.position = "top"))
ggsave("./plots/dfe_scatter_SE_v_SM_mutualism.png", plot = d,
       dpi = 300, width = 6, height = 6)

(e = fitness %>%
    mutate(treatment = ifelse(!E & !M, "monoculture",
                              ifelse(!E & M, "SM",
                                     ifelse(E & !M, "SE", "SEM")))) %>%
    group_by(locusId, treatment) %>%
    summarize(fitness = mean(fitness_normalized)) %>%
    pivot_wider(names_from = "treatment", values_from = "fitness") %>%
    left_join(results %>%
                mutate(locusId = gene) %>%
                select(term, locusId, padj) %>%
                pivot_wider(names_from = "term", values_from = "padj") %>%
                rename(intercept = `(Intercept)`)) %>%
    mutate(sign = ifelse(intercept < 0.05 & ETRUE >= 0.05 & MTRUE >= 0.05, "int.",
                         ifelse(intercept < 0.05 & ETRUE >= 0.05 & MTRUE < 0.05, "int. and main M",
                                ifelse(intercept < 0.05 & ETRUE < 0.05 & MTRUE >= 0.05, "int. and main E",
                                       ifelse(intercept < 0.05 & ETRUE < 0.05 & MTRUE < 0.05, "int. and both", "other"))))) %>%
    left_join(gene_names) %>%
    ggplot(aes(x = SE, y = SEM, color = sign))+
    scale_color_discrete(guide =  guide_legend(title = ""))+
    geom_abline(intercept = 0, slope = 1)+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_point()+
    ggrepel::geom_text_repel(size = 2.5, seed = 1, force = 0.5, force_pull = 1.5,
                             color = "black", segment.color = "gray",
                             aes(label = gene_name))+
    theme_bw(12)+
    theme(legend.position = "top"))
ggsave("./plots/dfe_scatter_SE_v_SEM_mutualism.png", plot = e,
       dpi = 300, width = 6, height = 6)

(f = fitness %>%
    mutate(treatment = ifelse(!E & !M, "monoculture",
                              ifelse(!E & M, "SM",
                                     ifelse(E & !M, "SE", "SEM")))) %>%
    group_by(locusId, treatment) %>%
    summarize(fitness = mean(fitness_normalized)) %>%
    pivot_wider(names_from = "treatment", values_from = "fitness") %>%
    left_join(results %>%
                mutate(locusId = gene) %>%
                select(term, locusId, padj) %>%
                pivot_wider(names_from = "term", values_from = "padj") %>%
                rename(intercept = `(Intercept)`)) %>%
    mutate(sign = ifelse(intercept < 0.05 & ETRUE >= 0.05 & MTRUE >= 0.05, "int.",
                         ifelse(intercept < 0.05 & ETRUE >= 0.05 & MTRUE < 0.05, "int. and main M",
                                ifelse(intercept < 0.05 & ETRUE < 0.05 & MTRUE >= 0.05, "int. and main E",
                                       ifelse(intercept < 0.05 & ETRUE < 0.05 & MTRUE < 0.05, "int. and both", "other"))))) %>%
    left_join(gene_names) %>%
    ggplot(aes(x = SM, y = SEM, color = sign))+
    scale_color_discrete(guide =  guide_legend(title = ""))+
    geom_abline(intercept = 0, slope = 1)+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_point()+
    ggrepel::geom_text_repel(size = 2.5, seed = 1, force = 0.5, force_pull = 1.5,
                             color = "black", segment.color = "gray",
                             aes(label = gene_name))+
    theme_bw(12)+
    theme(legend.position = "top"))
ggsave("./plots/dfe_scatter_SM_v_SEM_mutualism.png", plot = e,
       dpi = 300, width = 6, height = 6)
