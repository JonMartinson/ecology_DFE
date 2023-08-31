############
# This is the file containing functions for use with
# jmc_mult_t0_pipeline_20220722.r
#
# Many of these functions are re-writes of code from the Wetmore (et al. 2015)
# pipeline or the Morin (et al. 2019) pipeline.
#
# The main difference here is that we do sequence-depth correction before calculating
# fitness (like the Morin, but not Wetmore pipeline does) whereas otherwise we largely
# follow the Morin pipeline for filtering "bad" barcodes/genes/samples
#
############

library(tidyverse)

get_replicates_with_large_left_right_fitness_diff = function(fitness, threshold){
  # this returns a vector of replicate names if the median absolute difference
  # in fitness from the left-side and right-side of the gene is greater than the threshold
  fitness %>%
    group_by(locusId, replicate) %>%
    summarize(fitness_left = first(fitness_left),
              fitness_right = first(fitness_right)) %>%
    mutate(left_right_diff = abs(fitness_left - fitness_right)) %>%
    filter(!is.na(left_right_diff)) %>% # sometimes there arent reads on the left or right, not a lot, but some
    ungroup() %>%
    group_by(replicate) %>%
    summarize(mad12 = median(left_right_diff)) %>%
    filter(mad12 > threshold) %>%
    select(replicate) %>% unlist() %>% unname()
}

get_genes_with_few_reads_in_any_half = function(counts, T0s, threshold){
  # this returns a list of the GENES (aka locusId) where any of the T0 samples
  # have a count on the left half, or the right half, of the gene below the threshold
  counts %>%
    mutate(left_half_of_gene = f <= 0.5) %>%
    select(barcode, locusId, left_half_of_gene, T0s) %>%
    pivot_longer(-c(barcode, left_half_of_gene, locusId), names_to = "rep", values_to = "count") %>% 
    group_by(locusId, left_half_of_gene, rep) %>%
    summarize(reads_per_side = sum(count)) %>%
    ungroup() %>%
    group_by(locusId) %>%
    summarize(below_threshold = any(reads_per_side < threshold)) %>%
    filter(below_threshold) %>% 
    select(locusId) %>%
    unlist() %>% unname
}

get_genes_with_few_reads_in_chosen_samples = function(counts, T0s, threshold){
  # this returns a list of the GENES (aka locusId) where any of the T0 samples
  # have a count below the threshold
  counts %>%
    select(locusId, T0s) %>%
    pivot_longer(-locusId, names_to = "rep", values_to = "count") %>% 
    group_by(locusId, rep) %>%
    summarize(reads_per_gene = sum(count)) %>%
    ungroup() %>%
    group_by(locusId) %>%
    summarize(below_threshold = any(reads_per_gene < threshold)) %>%
    filter(below_threshold) %>% 
    select(locusId) %>%
    unlist() %>% unname
}

get_barcodes_with_any_below_threshold = function(counts, T0s, threshold){
  # this returns the list of barcodes where any of the T0 samples
  # have a count below the threshold
  counts %>%
    select(barcode, T0s) %>%
    pivot_longer(-barcode, names_to = "rep", values_to = "count") %>% 
    group_by(barcode) %>%
    summarize(below_threshold = any(count < threshold)) %>%
    filter(below_threshold) %>% 
    select(barcode) %>%
    unlist() %>% unname
}

get_barcodes_at_gene_extremes = function(counts, min_loc = 0.1, max_loc = 0.9){
  # examines the f value of the barcode to make sure it is >= min_loc and
  # <= max_loc, and returns the barcodes that fail this test
  counts %>%
    filter((f < min_loc) | (f > max_loc)) %>%
    select(barcode) %>%
    unlist() %>%
    unname()
}

get_replicates_with_low_median_gene_reads = function(counts, threshold){
  counts %>%
    pivot_longer(-base_columns, names_to = "replicate", values_to = "count") %>% 
    group_by(replicate, locusId) %>%
    summarize(reads_per_gene = sum(count)) %>%
    ungroup() %>%
    group_by(replicate) %>%
    summarize(median_reads_per_gene = median(reads_per_gene)) %>%
    ungroup() %>%
    filter(median_reads_per_gene < threshold) %>%
    select(replicate) %>%
    unname() %>% unlist()
}

add_pseudocount = function(counts, pseudocount = 0.1){
  counts %>%
    pivot_longer(cols = c(-barcode, -rcbarcode, -scaffold, 
                          -strand, -pos, -locusId, -f),
                 names_to = "replicate", values_to = "count") %>%
    mutate(raw_count = count, 
           count = count + pseudocount) 
}

filter_barcodes = function(counts, T0_thresh = 0, f_min = 0, f_max = 1){
  counts %>%
    filter(T0 > T0_thresh & f > f_min & f < f_max)
}

normalize_using_reference_genes = function(counts, ref){
  # this is NOT exaclty like the Dutton pipeline,
  # because it standardizes the normalized counts so that the
  # reference genes have the same median counts as the raw data
  # (rathre than having a median of 1)
  # this is important for using the counts to get variances later
  # based upon poisson assumptions
  counts %>%
    left_join(counts %>%
    filter(locusId %in% ref) %>%
    mutate(overall_mean_ref_count = mean(raw_count)) %>%
    group_by(replicate) %>%
    summarize(norm_factor = first(overall_mean_ref_count) / mean(raw_count),
              mean_ref_count = mean(raw_count),
              overall_mean_ref_count = first(overall_mean_ref_count))) %>%
    mutate(count_normalized = count * norm_factor)
}
  
calculate_strain_fitness_from_T0 = function(counts, T0_name = "T0"){
  # this calculates strain fitness and strain fitness variance
  # it uses the depth-normalized counts for fitness, and
  # it uses the raw counts for variance
  counts %>%
    mutate(log2count = log2(count_normalized)) %>%
    group_by(barcode) %>%
    mutate(strain_fitness = log2count - log2count[replicate == T0_name]) %>%
    mutate(strain_fitness_var = ((1 / (raw_count + 1)) + (1 / (raw_count[replicate == T0_name] + 1))) / log(2)^2) %>%
    ungroup() 
}

calculate_weighted_gene_fitness = function(fitness, count_cap = 20){
  # this calculates gene fitness by computing a weighted 
  # mean of the strain-level fitness, weighted by the inverse 
  # of their variance. Note, if the variance is above that which would
  # be found with a t1 and t0 each with 20 counts, the variance is capped
  # at that which would be found with a t1 and a t0 with 20 counts,
  # as in Wetmore et al. 
  #    This function also returns the weighted fitness on the left, and on the right
  #   halfves of the gene, for gene-wise QC
  cap = 1 / (((1/(count_cap+1))+(1/(count_cap+1)))/(log(2)^2))
  fitness %>% 
    mutate(strain_weight = 1/strain_fitness_var) %>%
    mutate(strain_weight_capped = ifelse(strain_weight > cap, 
                                         cap, strain_weight)) %>%
    group_by(replicate, locusId) %>%
    mutate(fitness = sum(strain_weight_capped*strain_fitness)/sum(strain_weight_capped),
              fitness_left = sum(strain_weight_capped[f <= 0.5] *strain_fitness[f <= 0.5])/sum(strain_weight_capped[f <= 0.5]),
              fitness_right = sum(strain_weight_capped[f > 0.5] *strain_fitness[f > 0.5])/sum(strain_weight_capped[f > 0.5])) %>%
    ungroup() 
}

calculate_fitness_from_T0 = function(counts){
  counts %>%
    pivot_longer(cols = c(-barcode, -rcbarcode, -scaffold, 
                          -strand, -pos, -locusId, -f),
                 names_to = "replicate", values_to = "count") %>%
    mutate(log2count = log2(count)) %>%
    group_by(barcode, locusId) %>%
    mutate(strain_fitness = log2count - log2count[replicate == "T0"]) %>%
    ungroup() %>%
    filter(replicate != "T0") %>%
    group_by(replicate, locusId) %>%
    summarize(fitness_var = var(fitness),
              fitness = mean(fitness)) %>%
    ungroup() 
}


subtract_rolling_median = function(fitness, genes.tab, window_size = 251){
  fitness_corrected = data.frame()
  scaffolds = unique(fitness$scaffold)
  gene_fitness = left_join(fitness %>%
                        group_by(replicate, locusId, scaffold) %>%
                        summarize(fitness = first(fitness)), 
                      genes.tab, by = c("locusId")) 
  for (curr_scaffold in scaffolds){
    n_genes = length(unique(gene_fitness$locusId[gene_fitness$scaffold == curr_scaffold]))
    if (n_genes >= window_size){
      this_scaffold = gene_fitness %>%
        filter(scaffold == curr_scaffold) %>%
        group_by(replicate) %>%
        arrange(begin) %>%
        mutate(fitness_normalized = fitness - runmed(fitness, window_size,endrule = "constant"))
    }else if (n_genes > 10){
      this_scaffold = gene_fitness %>%
        filter(scaffold == curr_scaffold) %>%
        group_by(replicate) %>%
        arrange(begin) %>%
        mutate(fitness_normalized = fitness - median(fitness))    
    }else{
      this_scaffold = gene_fitness %>%
        filter(scaffold == curr_scaffold) %>%
        group_by(replicate) %>%
        arrange(begin) %>%
        mutate(fitness_normalized = fitness)        
    }
    fitness_corrected = rbind(fitness_corrected, 
                              data.frame(this_scaffold)) %>%
      select(-begin, -fitness)
  }
  left_join(fitness, fitness_corrected, by = c("locusId",
                                               "replicate")) 
}


############ Novel functions (not just rewrites of Dutton code)
analyze_fitness = function(fitness, response_var = "fitness_gen_corr"){
  # does the two-factor lm on each gene and saves the data, then does BH corrections
  # it does this on 
  results = data.frame()
  for (gene in unique(fitness$locusId)){
    fitness %>%
      filter(locusId == gene) %>%
      rename_(response_var = "response_var") %>%
      lm(response_var ~ S + M + S:M, data = .) %>%
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



