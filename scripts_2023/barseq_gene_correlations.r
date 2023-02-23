rm(list = ls())
library(tidyverse)
library(ggridges)
options(dplyr.summarise.inform = FALSE)

# get_similar = function(fitness_matrix, base_gene, thresh){
#   similar = names(fitness_matrix[base_gene,])[abs(fitness_matrix[base_gene,]) > thresh]
#   cors = fitness_matrix[base_gene,][abs(fitness_matrix[base_gene,]) > thresh]
#   data.frame(locusId = similar, correlation = cors)
# }
# # different attempt at cofitness. One prob with pearson correlation is that
# # it doesn't care about magnitude, or have any info on differences between treatments
# # Here I add a few datapoints--the distances between each treatment
# fitness_summary = fitness %>%
#   filter(ecology == "mutualism") %>%
#   group_by(community,locusId) %>%
#   summarize(fitness = mean(fitness_normalized)) %>%
#   ungroup() %>%
#   mutate(estimate = community) %>%
#   select(estimate, locusId, fitness) %>%
#   pivot_wider(names_from = locusId, values_from = fitness) 
# 
# signed_dists = function(x){
#   dists = rep(0, (length(x) *(  length(x)-1))/ 2)
#   counter = 1
#   for (i in 1:(length(x)-1)){
#     for (j in (i+1):length(x)){
#       dists[counter] = x[i] - x[j]
#       counter = counter + 1
#     }
#   }
#   dists
# }
# dists = apply(fitness_summary[,2:ncol(fitness_summary)], 2, 
#               FUN = signed_dists) %>%
#   as.data.frame() %>%
#   mutate(estimate = paste("dist", 1:6))
# fitness_summary = rbind(fitness_summary, dists)
# fitness_matrix = matrix(0, nrow = length(unique(fitness$locusId)),
#                         ncol = length(unique(fitness$locusId)),
#                         dimnames = list(unique(fitness$locusId),
#                                         unique(fitness$locusId)))
# for (gene1 in unique(fitness$locusId)){
#   print(gene1)
#   for (gene2 in unique(fitness$locusId)){
#     fitness_matrix[gene1, gene2] = cor(fitness_summary[gene1], fitness_summary[gene2])
#   }
# }
# fitness_matrix_dists = fitness_matrix
# save(fitness_matrix_dists, file = "./data/S0_barseq_wetlab/fitness_correlations_mutualism_with_dists.Rdat")

# calculates "cofitness" (see Price et al 2018) and has analyses

fitness = read_csv("./data/S0_barseq_fitnesses_20220722.csv")
gene_names = read_csv("./data/uniprot_s_enterica_gene_names.csv") %>%
  mutate(multi_name = `Gene Names`,
         gene_name = `Gene Names (primary)`) %>%
  mutate(locusId = sapply(multi_name,
                           FUN = function(x){
                             str_split(x, " ")[[1]][length(str_split(x, " ")[[1]])]
                           })) %>%
  select(gene_name, locusId)
#
mutualism <- c("T2 S mut", "T3 SE mut", "T5 SM mut", "T5 SEM mut") # early mutualism samples
competition = c("T6 S comp", "T6 SE comp", "T6 SM comp", "T6 SEM comp")
all_trts = c(mutualism, competition)
#
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
                                   ifelse(M, "SM", "S")))) %>%
  mutate(community = factor(community,
                            levels = c("S","SE","SM","SEM")))
#
fitness_summary = fitness %>%
  group_by(community, ecology, locusId) %>%
  summarize(fitness = mean(fitness_normalized)) %>%
  ungroup() %>%
  select(community, ecology, locusId, fitness) %>%
  pivot_wider(names_from = locusId, values_from = fitness)

fitness_matrix = matrix(0, nrow = length(unique(fitness$locusId)),
                     ncol = length(unique(fitness$locusId)),
                     dimnames = list(unique(fitness$locusId),
                                     unique(fitness$locusId)))
for (gene1 in unique(fitness$locusId)){
  print(gene1)
  for (gene2 in unique(fitness$locusId)){
 fitness_matrix[gene1, gene2] = cor(fitness_summary[gene1], fitness_summary[gene2])
  }
}
save(fitness_matrix, file = "./data/fitness_correlations.Rdat")

# lets repeat it, but also add a point per locusId which is the mutualism variance,
# and a point which is the competition variance.
#
# fitness_summary = fitness %>%
#   group_by(community, ecology, locusId) %>%
#   summarize(fitness = mean(fitness_normalized)) %>%
#   ungroup() %>%
#   select(community, ecology, locusId, fitness) %>%
#   bind_rows(fitness %>%
#               group_by(community, ecology, locusId) %>%
#               summarize(fitness = mean(fitness_normalized)) %>%
#               ungroup() %>%
#               group_by(ecology, locusId) %>%
#               summarize(community = "sd", fitness = sd(fitness))) %>%
#   pivot_wider(names_from = locusId, values_from = fitness)
# #
# fitness_matrix_sd = matrix(0, nrow = length(unique(fitness$locusId)),
#                         ncol = length(unique(fitness$locusId)),
#                         dimnames = list(unique(fitness$locusId),
#                                         unique(fitness$locusId)))
# for (gene1 in unique(fitness$locusId)){
#   print(gene1)
#   for (gene2 in unique(fitness$locusId)){
#     fitness_matrix_sd[gene1, gene2] = cor(fitness_summary[gene1], fitness_summary[gene2])
#   }
# }
#lets repeat it, but only correlate the mutualism data
#
# fitness_summary = fitness %>%
#   group_by(community, ecology, locusId) %>%
#   summarize(fitness = mean(fitness_normalized)) %>%
#   ungroup() %>%
#   filter(ecology == "mutualism") %>%
#   select(community, locusId, fitness) %>%
#   pivot_wider(names_from = locusId, values_from = fitness)
#
# fitness_matrix_m = matrix(0, nrow = length(unique(fitness$locusId)),
#                         ncol = length(unique(fitness$locusId)),
#                         dimnames = list(unique(fitness$locusId),
#                                         unique(fitness$locusId)))
# for (gene1 in unique(fitness$locusId)){
#   print(gene1)
#   for (gene2 in unique(fitness$locusId)){
#     fitness_matrix_m[gene1, gene2] = cor(fitness_summary[gene1], fitness_summary[gene2])
#   }
# }
# #save(fitness_matrix_m, file = "./data/S0_barseq_wetlab/fitness_correlations_only_mutualism.Rdat")
# save(fitness_matrix_sd, file = "./data/S0_barseq_wetlab/fitness_correlations_with_sd.Rdat")
save(fitness, file = "./data/fitness.Rdat")
load("./data/S0_barseq_wetlab/fitness_correlations_only_mutualism.Rdat")                        
load("./data/S0_barseq_wetlab/fitness_correlations_with_sd.Rdat")                        
load("./data/S0_barseq_wetlab/fitness_correlations.Rdat")                        
load("./data/S0_barseq_wetlab/fitness.Rdat")
load("./data/S0_barseq_wetlab/genenames.Rdat")
results = read_csv("./data/S0_barseq_wetlab/mutualism_two_way_anova_stats.csv")

gene_names$locusId[gene_names$gene_name == "ilvA"] = "STM3905"
gene_names$locusId[gene_names$gene_name == "ilvC"] = "STM3909"
gene_names$locusId[gene_names$gene_name == "ilvD"] = "STM3904"
gene_names$locusId[gene_names$gene_name == "ilvE"] = "STM3903"

gene_names$gene_name[gene_names$locusId == "STM3531"] = "STM3531"
gene_names$gene_name[gene_names$locusId == "STM2389"] = "fadI"
gene_names$gene_name[gene_names$locusId == "STM3982"] = "fadA"
gene_names$gene_name[gene_names$locusId == "STM3983"] = "fadB"

save(gene_names, file = "./data/gene_names.Rdat")


# overall distribution of correlations
hist(fitness_matrix, 20)

# correlations of gene STM0774 with others
base_gene = "STM2501" # ppk
base_gene = "STM2338" # pta
base_gene = "STM2337" # ackA
base_gene = "STM4275" # acs
base_gene = "STM0182" # panB
base_gene = "STM3901" # ilvG https://ecocyc.org/ECOLI/NEW-IMAGE?type=PATHWAY&object=BRANCHED-CHAIN-AA-SYN-PWY
base_gene = "STM0959" # lrp, which represses ilvG/M
base_gene = "STM1339" # ihfA, which activates ilvG/M. this is lowest fitness in SEM
base_gene = "STM0113" # leuA, which is downstream of ilvG/M and is in leucine biosynthesis
base_gene = "STM0112" # leuB 
base_gene = "STM1883" # related to ackA, makes ATP from acetyl phosphate to make acetate (or reverse)
base_gene = "STM3905" # ilvA
base_gene = "STM4275" # acs
base_gene = "STM3796" # ilvB
base_gene = "STM3909" # ilvC
base_gene = "STM0736" # sucA
base_gene = "STM1818" # fadD
base_gene = "STM0309" # fadE

hist(fitness_matrix[base_gene,], 20)



gene_names %>%
  filter(str_detect(gene_name, "ilv"))

#similar genes by correlation threshold
thresh = 0.98

similar_no_sd = get_similar(fitness_matrix_dists, base_gene, thresh)

fitness %>%
  filter(locusId %in% similar_no_sd$locusId) %>%
  left_join(gene_names) %>%
  left_join(similar_no_sd) %>%
  filter(correlation > 0) %>%
  filter(ecology == "mutualism") %>%
  mutate(locus = paste(gene_name, locusId)) %>%
  mutate(locus = factor(locus,
                          levels = unique(locus[order(correlation)]))) %>%
  ggplot(aes(x = community, y = fitness_normalized))+
  geom_point(shape = 1)+
  stat_summary(fun = mean, shape = "-", size = 2)+
  facet_wrap(~locus, scales = "free_y")+
  theme_bw(8)

similar_sd = get_similar(fitness_matrix_sd, base_gene, thresh = 0.9)

fitness %>%
  filter(locusId %in% similar_sd$locusId) %>%
  left_join(gene_names) %>%
  left_join(similar_sd) %>%
  filter(correlation > 0) %>%
  filter(ecology == "mutualism") %>%
  mutate(locusId = factor(locusId,
                          levels = unique(locusId[order(correlation)]))) %>%
  ggplot(aes(x = community, y = fitness_normalized))+
  geom_point(shape = 1)+
  stat_summary(fun = mean, shape = "-", size = 2)+
  facet_wrap(gene_name~locusId, nrow = 2)

fitness %>%
  left_join(gene_names) %>%
  #filter(str_detect(gene_name, "ilv") | 
  #         str_detect(gene_name, "pkjkgm")) %>%
  #filter(gene_name == "gltA") %>%
  filter(locusId == "STM1378") %>%
  filter(ecology == "mutualism") %>%
  mutate(locus = paste(gene_name, locusId)) %>%
  ggplot(aes(x = community, y = fitness_normalized))+
  geom_point(shape = 1)+
  stat_summary(fun = mean, shape = "-", size = 2)+
  facet_wrap(~locus, scales = "free_y")+
  theme_bw(8)

# genes from FBA which can split needed during growth in pyruvate vs. maltose
pyr = c('STM0999', 'STM0320', 's0001', 'STM0736', 'STM0154',
        'STM0733', 'STM4109', 'STM4290', 'STM0734', 'STM1473', 
        'STM0517', 'STM1349', 'STM0732', 'STM2267', 'STM0735', 'STM0737')
malt = c('STM4115', 'STM0999', 'STM0694', 'STM1749', 'STM1480',
         'STM3957', 'STM4227', 'STM2461', 'STM3661', 'STM3514', 
         'STM1651', 'STM4341', 'STM3692', 'STM4340', 'STM0320',
         'STM4231', 'STM0488', 'STM1888', 'STM1710', 'STM0984',
         'STM4229', 'STM3045', 'STM2526', 'STM4084', 'STM4086',
         'STM1479', 'STM2463', 'STM1567', 'STM1378', 'STM1627',
         'STM4568', 'STM4062', 'STM4230', 'STM0970', 'STM2267',
         'STM1326', 'STM4228', 'STM4451', 'STM0973', 'STM0420',
         'STM1700', 'STM4343', 'STM4570', 'STM1473', 'STM4114',
         'STM0844', 'STM3241', 'STM4342', 'STM4452', 'STM0843')

fitness %>%
  left_join(gene_names) %>%
  filter(locusId %in% malt) %>%
  filter(ecology == "mutualism") %>%
  mutate(locus = paste(gene_name, locusId)) %>%
  ggplot(aes(x = community, y = fitness_normalized))+
  geom_point(shape = 1)+
  stat_summary(fun = mean, shape = "-", size = 2)+
  facet_wrap(~locus, scales = "free_y")+
  theme_bw(8)

# genes which are low fitness in presence of E
E_low = results %>%
  filter(term == "ETRUE") %>%
  filter(padj < 0.05 & Estimate < 0) %>%
  select(gene) %>% pull 

# genes which are improved only in presence of E
E_low = results %>%
  group_by(gene) %>%
  filter(padj[term == ""])
  filter(term == "ETRUE") %>%
  filter(padj < 0.05 & Estimate < 0) %>%
  select(gene) %>% pull 

fitness %>%
  left_join(gene_names) %>%
  filter(locusId %in% E_low) %>%
  filter(ecology == "mutualism") %>%
  mutate(locus = paste(gene_name, locusId)) %>%
  ggplot(aes(x = community, y = fitness_normalized))+
  geom_point(shape = 1)+
  stat_summary(fun = mean, shape = "-", size = 2)+
  facet_wrap(~locus, scales = "free_y")+
  theme_bw(8)

fitness %>%
  left_join(gene_names) %>%
  filter(gene_name %in% c("envZ", "ompF", "ppk")) %>%
  filter(ecology == "mutualism") %>%
  mutate(locus = paste(gene_name, locusId)) %>%
  ggplot(aes(x = community, y = fitness_normalized))+
  geom_point(shape = 1)+
  stat_summary(fun = mean, shape = "-", size = 2)+
  facet_wrap(~locus, scales = "free_y")+
  theme_bw(8)+
  labs(x = "", y = "fitness")
Cfind_correlation_family = function(fitness_matrix_sd, base_gene, thresh = 0.95){
  # iteratively gets all the genes with correlation (or anticorrelation)
  # higher than thresh.
  family = c(base_gene)
  new_members = names(fitness_matrix_sd[base_gene,])[abs(fitness_matrix_sd[base_gene,]) > thresh]
  while (length(setdiff(new_members, family)) > 0){
    family = c(family, new_members)
    newest_members = c()
    for (member in new_members){
      newest_members = c(newest_members,
                         names(fitness_matrix_sd[member,])[abs(fitness_matrix_sd[member,]) > thresh])
      
    }
    new_members = unique(newest_members)
  }
  return(unique(family))
}
family = find_correlation_family(fitness_matrix_sd, base_gene)


fitness %>%
  filter(ecology == "mutualism") %>%
  filter(locusId %in% new_members) %>%
  ggplot(aes(x = community, y = fitness_normalized))+
  geom_point()+
  stat_summary(fun = mean, shape = "-", size = 2)+
  facet_wrap(~locusId)
