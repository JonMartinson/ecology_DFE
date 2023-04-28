###################
#
# This is the current pipeline with all filtering,
# up to the point where fitnesses are measured. 
#
###################


rm(list = ls())
library(tidyverse)

# file location formatting may change depending on the OS you are running. 
# Current locations format: MacOS
source("scripts_2023/pipeline/functions/pipeline_functions.r") # location of pipeline functions
scaffoldX=  "AE006468" # Salmonella LT2 chromosome 
ref=c("STM0604","STM1237","STM0329","STM2774","STM0333") # here you write the neutral gene references



gene_tab_loc = "data/genes.GC.txt" #GC content of genes, contains gene location etc.
harcpoolcount_loc = "data/all.poolcount" #output from feba barseq pipeline Wetmore 2015
metadata_loc = "data/So_5_23_22_PCR_layout.csv" # identities of barcodes, etc.


# Dutton parameter values were 0.1, 3.1, 0.1, 0.9, 251
pseudocount = 0.1 # what to add before doing log2() ?
T0_thresh = 3.1  # minimum count for a T0 barcode
f_min = 0.1 # minimum position for a barcode
f_max = 0.9 # maximum position for a barcode
window_size = 251 # num genes to subtract a smoothed median from 
# this is needed to go through the workflow per-replicate
T0s = paste("so_apr_22.IT00", 1:5, sep = "")


# note:  the data in the file in harcpoolcount_loc must have columns called "barcode", "rcbarcode", "scaffold",
# "strand", "pos","locusId",  "f", "T0" (and then 1 or more treatment columns)
base_columns = c("barcode","rcbarcode","scaffold","strand","pos","locusId",
                 "f")


#####read files

metadata = read_csv(metadata_loc) %>%
  mutate(replicate = paste("so_apr_22", `Primer P2`, sep = "."),
         treatment = `Sample Identity`) %>%
  filter(!is.na(treatment)) %>% # because the excel sheet had 5 extra rows
  select(replicate, treatment) %>%
  mutate(treatment = ifelse(treatment == "T0 S inoculum", "T0", treatment)) %>%
  group_by(treatment) %>%
  mutate(rep = 1:n())

# gene data
genes.tab <- readr::read_delim(gene_tab_loc, "\t", escape_double = FALSE, trim_ws = TRUE)
# all count data--no need to do this piecemeal like in Dutton's scripts
raw_counts = read.csv(harcpoolcount_loc, sep = "\t", stringsAsFactors = FALSE) 


#########
# Data filtering
######## 

# 1. Only keep reads in the chromosome scaffold
counts = raw_counts %>%
  filter(scaffold == scaffoldX) 
print(paste("removed", 
            100 * (1 - (nrow(counts) / nrow(raw_counts))),
            " percent of barcodes not on the chromosome"))

# 2. remove the chromosomal barcodes which are not in a gene
print(paste("removed", 
            100 * ((sum(counts$locusId == "") / nrow(counts))),
            " percent of remaining barcodes because they are not in a gene"))
counts = counts %>%
  filter(locusId != "") 

# 3. remove barcodes in the first 10% or last 10% of a gene
extreme_loc_barcodes = get_barcodes_at_gene_extremes(counts, f_min, f_max)
print(paste("removing", 
            100 * ((length(extreme_loc_barcodes) / nrow(counts))),
            " percent of remaining barcodes because they are at the extreme ends of a gene"))
counts = counts %>%
  filter(!(barcode %in% extreme_loc_barcodes))


# 4. remove barcodes which never have more than 3 reads in the T0 samples
low_barcodes = get_barcodes_with_any_below_threshold(counts, T0s, 4) # this is what Dutton used (text says 3, but code says 4 due to >)
print(paste("removed", 
            100 * ((length(low_barcodes) / nrow(counts))),
            " percent of remaining barcodes because the read count is too low in at least one T0"))
counts = counts %>%
  filter(!(barcode %in% low_barcodes))

# 5. Identify samples with less than a median of 50 reads per gene 
low_gene_count_reps = get_replicates_with_low_median_gene_reads(counts, 50)
print(paste("removing reps: ", low_gene_count_reps, ", due to low median gene counts [if empty, all passed threshold]"))
counts = counts %>%
  select(-low_gene_count_reps)

# 6. Remove GENES which have less than 30 reads across all barcodes
low_genes = get_genes_with_few_reads_in_chosen_samples(counts, T0s, 30)
print(paste("removed", 
            100 * ((length(low_genes) / length(unique(counts$locusId)))),
            " percent of remaining genes because the gene-wise read count is too low in the T0s"))
counts = counts %>%
  filter(!(locusId %in% low_genes))

# 7. remove genes with less than 15 reads on left, or right half of gene
low_gene_halves = get_genes_with_few_reads_in_any_half(counts, T0s, 15)
print(paste("removed", 
            100 * ((length(low_gene_halves) / length(unique(counts$locusId)))),
            " percent of remaining genes because the left-half or right-half gene read count is too low in the T0s"))
counts = counts %>%
  filter(!(locusId %in% low_gene_halves))

print(paste("total barcodes removed = ", nrow(raw_counts) - nrow(counts),
            " (", (nrow(raw_counts) - nrow(counts)) / nrow(raw_counts), "%)"))

length(unique(counts$locusId))
length(unique(genes.tab$locusId))
length(unique(raw_counts$locusId))
### the main loop; it uses a different T0 for each strain fitness calculation,
### doing it per-rep. An alternative would be to use an average T0 as the base
### for each fitness calculation.

strain_fitness = data.frame()
for (rep in unique(metadata$rep)){
  print(rep)
  # subsample this rep's data, with a different T0 per loop
  samples = metadata$replicate[metadata$rep == rep]
  T0_sample = metadata$replicate[metadata$rep == rep & metadata$treatment == "T0"]
  this_rep = counts %>%
    select(c(base_columns, samples)) %>%
    rename(T0 = T0_sample)
  
  # perform the rest of the fitness calculation for this replicate
  this_rep = add_pseudocount(this_rep, 0.1)
  this_rep = normalize_using_reference_genes(this_rep, ref)
  this_rep = calculate_strain_fitness_from_T0(this_rep) %>%
    mutate(replicate = ifelse(replicate == "T0", T0_sample, replicate))

  strain_fitness = rbind(strain_fitness, this_rep)
  rm(this_rep)
}
fitness = strain_fitness %>%
  filter(!(replicate %in% T0s))

# count cap makes it so that barcodes with very high read counts (greater than count_cap)
# don't have outsized leverage on the mean
fitness = calculate_weighted_gene_fitness(fitness, count_cap = 20)

n_samples = fitness %>%
  left_join(metadata) %>%
  filter(treatment %in% c("T2 S mut", "T3 SE mut", "T5 SM mut", "T5 SEM mut")) %>%
  select(replicate) %>% unlist() %>% unname() %>% unique() %>% length()

# Find genes which have median fitness on left / right that is too uneven
too_uneven_reps = get_replicates_with_large_left_right_fitness_diff(fitness, 0.5)
print(paste("removing reps: ", paste(too_uneven_reps, collapse = ", "), ", due to uneven left/right fitness [if empty, all passed threshold]"))
fitness = fitness %>%
  filter(!(replicate %in% too_uneven_reps))


### remove chromosome location bias
fitness = subtract_rolling_median(fitness, genes.tab)
rm("counts","strain_fitness","raw_counts", "extreme_loc_barcodes",
   "low_barcodes", "low_gene_count_reps","low_gene_halves","low_genes")
##
# At this point, the location-normalized fitnesses are calculated ("fitness_normalized")
# IMPORTANT: all barcode-level data is still retained, so when one does 
# stats on the fitnesses, one must first summarize down to a single value 
# per replicate-locusId (as seen below)
##

fitness %>%
  left_join(metadata) %>%
  group_by(replicate, treatment, locusId) %>%
  summarize(fitness_normalized = first(fitness_normalized)) %>%
  ggplot(aes(x = fitness_normalized, fill = treatment))+
  geom_density()+
  geom_vline(xintercept = 0, linetype = "dashed")+
  facet_wrap(~replicate)

write_csv(fitness %>% left_join(metadata),
          "data/S0_barseq_fitnesses_20220722.csv")


