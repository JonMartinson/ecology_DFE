rm(list = ls())

library(tidyverse)
library(lubridate)
library(rstatix)

#### initial data wrangling ####

df_count <- read_csv("data/sample_counts.csv") # plate count data
df_OD <-
    read_csv("data/time_density.csv") # OD data and dates/times that sampling occurred

df <- left_join(df_count, df_OD)

## mutualistic S monoculture rep # 1 did not grow so I eliminated it from the data frame and
## replaced it with mutualistic S monoculture rep # 6

df <-  df %>% 
    filter(!(sample_names == 'T2 S mut' & is.na(glu_count))) %>% 
    mutate(rep = case_when(
        sample_names == 'T2 S mut' & rep == 6 ~ 1,
        TRUE ~ rep))

# day and hour in which the experiment started
time_zero <-
    ymd_hms("2022-04-14 19:00:00")

# mL cells plated (5 µL drops) for enumeration (plate counts)
amt_plated <- 0.005 

# function to calculate cfu/mL for each species
CFU_mL <- function(count, dil) {
        count / (dil * amt_plated)
    }



df <- df %>%
    separate(
        col = sample_names,
        into = c("sampling_timepoint", "species_present", "ecology"),
        remove = FALSE
    ) %>%
    unite("date_time_0", c(Day, Hour), sep = " ") %>% # combine date and time columns
    mutate(date_time = mdy_hms(date_time_0)) %>%  # make date_time_0 parsable with lubridate
    select(-date_time_0) %>% # remove temporary column
    mutate(hours_incubated = date_time - time_zero) %>% # calculate amt of time from T0 to sampling time
    mutate(
        S_CFU_mL = CFU_mL(glu_count, glu_dil),
        # calculate CFU/mL for each species
        E_CFU_mL = CFU_mL(lac_met_count, lac_met_dil),
        M_CFU_mL = CFU_mL(succ_count, succ_dil)
    )

# get the population sizes for each species at T0
t0_pop_size <- df %>% 
    mutate(
        total_mL_onto_plate = # total mL inoculated onto experimental plates for each species combination
            case_when(
                species_present == 'S' ~ 0.025, # 0.025mL or 25µL pipetted onto experimental plate
                species_present == 'SE' ~ 0.05,
                species_present == 'SM' ~ 0.05,
                species_present == 'SEM' ~ 0.075
            )
    ) %>%
    filter(sampling_timepoint == 'T0') %>%
    group_by(species_present) %>%
    mutate(
        S_T0_pop_size = S_CFU_mL * total_mL_onto_plate, # S Cfu/mL * mL inoculated onto experimental plate
        E_T0_pop_size = E_CFU_mL * total_mL_onto_plate,
        M_T0_pop_size = M_CFU_mL * total_mL_onto_plate
    ) %>%
    summarise(
        T0_S_pop = mean(S_T0_pop_size),
        T0_E_pop = mean(E_T0_pop_size),
        T0_M_pop = mean(M_T0_pop_size)
    ) 



## add starting population sizes for each species for each species combination
## then calculate the number of generations for each species for each experimental replicate
## the CFU/mL is multiplied by 4 because 4 mL of saline was used to scrape the cells off of the plate
df <- df %>% 
    left_join(t0_pop_size, by = 'species_present') %>% 
    mutate(S_generations = log2((S_CFU_mL*4)/T0_S_pop), # multiplied CFU/mL * 4 because 4mL of cells were scraped up 
           E_generations = log2((E_CFU_mL*4)/T0_E_pop),
           M_generations = log2((M_CFU_mL*4)/T0_M_pop)) %>% 
    mutate(ecology = case_when(ecology == 'mut' ~ 'Mutualism', # change the names of ecology to Mutualism and competition
                               ecology == 'comp' ~ 'Competition'
    ))


#### Generations ####
## number of generations with both mutualism and competitive interactions
## line goes through the middle of the mean of each treatment timepoint

df %>% 
    filter(str_detect(ecology,"Mutualism|Competition")) %>% 
    group_by(species_present, ecology) %>% 
    group_modify(~ add_row(.x,hours_incubated = duration(num = 0, units = "seconds"), S_generations = 0)) %>% 
    ggplot(aes((hours_incubated/60)/60, S_generations, color = factor(species_present, levels = c('S', 'SE', 'SM', 'SEM')))) +
    stat_summary(fun.y=mean, geom="line")+
    geom_jitter()+
    viridis::scale_color_viridis(discrete = T) +
    labs(title = "Number of Salmonella Generations") +
    xlab("Hours of Incubation")+
    ylab("# of Generations") +
    theme_bw(12) +
    labs(color = 'Community')+
    facet_wrap(~factor(ecology, levels = c('Mutualism', 'Competition')))

ggsave("plots/number_of_S_gens_through_time.png",
       dpi = 300, width = 8, height = 4)

## look just at the number generations in the treatments that we performed barseq on

target <- c("T2 S mut", "T3 SE mut","T5 SM mut","T5 SEM mut", 'S comp',"E comp", "M comp") %>% 
    str_c( collapse = "|")



generations_species <- df %>% 
    filter(str_detect(sample_names, target)) %>% 
    select(sample_names, sampling_timepoint,rep, ecology ,species_present ,S_generations, E_generations, M_generations) %>% 
    filter(rep != 6) %>% 
    group_by(sample_names) %>% 
    pivot_longer(cols = c(S_generations, E_generations, M_generations), 
                 names_to = "species",
                 values_to = "generations") %>% 
    ungroup() %>% 
    mutate(species = case_when(
        species == 'S_generations' ~ 'S',
        species == 'E_generations' ~ 'E',
        species == 'M_generations' ~ 'M')) %>% 
    ggplot(aes(x = factor(species_present, levels = c('S','SE', 'SM', 'SEM')), y = generations))+
    stat_summary(fun = mean, shape = '-', size = 3, color = 'black')+
    geom_point(shape = 1)+
    theme_bw(12)+
    xlab("Treatment")+
    ylab("Generations")+
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = 'none') +
    facet_grid(factor(species, levels = c('S','E', 'M'))~ factor(ecology, levels = c('Mutualism', 'Competition')))

saveRDS(generations_species, file = 'rds_plots/generations_species.rdata')

ggsave("plots/number_of_species_generationss_for_experimental_treatments.png",
       dpi = 300, width = 8, height = 4)

#### frequency of species ####

df_freq <- df %>% 
    mutate(
        S_CFU_mL = replace_na(S_CFU_mL, 0),
        E_CFU_mL = replace_na(E_CFU_mL, 0),
        M_CFU_mL = replace_na(M_CFU_mL, 0),
        freq_S = S_CFU_mL/(S_CFU_mL + E_CFU_mL + M_CFU_mL),
        freq_E = E_CFU_mL/(S_CFU_mL + E_CFU_mL + M_CFU_mL),
        freq_M = M_CFU_mL/(S_CFU_mL + E_CFU_mL + M_CFU_mL),
        S_total_CFU = S_CFU_mL*4, # *4 because 4 mL was used to scrape plates
        E_total_CFU = E_CFU_mL*4,
        M_total_CFU = M_CFU_mL*4
    ) 

#make long frequency data frame for plotting purposes
df_freq_long <- df_freq %>% 
    select(sample_names, sampling_timepoint,rep, ecology ,species_present ,freq_S, freq_E, freq_M) %>% 
    group_by(sample_names) %>% 
    pivot_longer(cols = c(freq_S, freq_E, freq_M), 
                 names_to = "species",
                 values_to = "species_frequency")

# make a filtered data frame showing just the treatments we performed barseq on

target <-
    c("T3 SE mut", "T5 SM mut", "T5 SEM mut", "E comp", "M comp") %>%
    str_c(collapse = "|")

df_freq_long_target <- df_freq_long %>% 
    filter(str_detect(sample_names, target))

df_freq_long_target$sample_names<- str_remove(df_freq_long_target$sample_names, 
                                              pattern = "^.{0,3}")


## show frequency of species in each replicate for all treatments

species_frequencies_end <- df_freq_long_target %>% 
    mutate(
           species = case_when(species == 'freq_E' ~ 'E',
                               species == 'freq_M' ~ 'M',
                               species == 'freq_S' ~ 'S')) %>% 
    filter(rep != 6) %>% #remove 6th replicate to avoid confusion since we only sequenced 5 reps
    ggplot(aes(x = rep, y = species_frequency, fill = factor(species, levels = c('S', 'E', 'M')))) +
    geom_col() + 
    labs(title = "Frequency of Species Engaged in Mutualism and Competition")+
    xlab("Replicate") +
    ylab("Species Frequency") +
    facet_grid(factor(species_present, levels = c('SE', 'SM', 'SEM'))~factor(ecology, levels = c('Mutualism', 'Competition')))+
    theme_bw(12)

saveRDS(species_frequencies_end, file = 'rds_plots/species_frequencies_end.rdata')

ggsave("plots/species_freqencies_all_reps.png",
       dpi = 300, width = 8, height = 4)

#### fitness of S in each condition ####

S_t0_freqs <-  df_freq_long %>% 
    filter(sampling_timepoint == 'T0') %>% 
    filter(species == 'freq_S') %>% 
    group_by(species_present) %>% 
    summarise(T0_S_freq_mean = mean(species_frequency)) %>% 
    ungroup()

S_fitness <- df_freq_long_target %>%
    left_join(S_t0_freqs, by = 'species_present') %>%
    filter(species == 'freq_S') %>%
    group_by(sampling_timepoint, rep, species_present, ecology) %>%
    summarise(fitness = log2(species_frequency) - log2(T0_S_freq_mean))

salmonella_fitness_barseq <- S_fitness %>% 
    filter(rep != 6) %>% 
    ggplot(aes(x = factor(species_present, levels = c('SE', 'SM', 'SEM')), y = fitness))+
    stat_summary(fun = mean, shape = '-', size = 3, color = 'black')+
    geom_point(shape = 1)+
    theme_bw(12)+
    xlab("Treatment")+
    ylab("Fitness")+
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = 'none') +
    facet_wrap(~ factor(ecology, levels = c('Mutualism', 'Competition'))) +
    ggtitle("Salmonella fitness in communities of differing ecology")

saveRDS(object = salmonella_fitness_barseq, file = 'rds_plots/salmonella_fitness_barseq.rdata' )

ggsave("plots/Salmonella_species_fitness_in_comp_mut_treatments.png",
       dpi = 300, width = 8, height = 4)

#### Total population size ####

#make long population size data frame for plotting purposes
df_total_long <- df_freq %>% 
    select(sample_names, sampling_timepoint,rep, ecology ,species_present ,S_total_CFU, E_total_CFU, M_total_CFU) %>% 
    group_by(sample_names) %>% 
    pivot_longer(cols = c(S_total_CFU, M_total_CFU, E_total_CFU), 
                 names_to = 'species',
                 values_to = 'species_total_cfu')

rm(target)
target <- c("T2 S mut", "T3 SE mut","T5 SM mut","T5 SEM mut", 'S comp',"E comp", "M comp") %>% 
    str_c( collapse = "|")

df_total_long_target <- df_total_long %>% 
    filter(str_detect(sample_names, target))

## all species
total_pop_size_barseq <- df_total_long_target %>%
    mutate(
           species = case_when(species == 'E_total_CFU' ~ 'E',
                               species == 'M_total_CFU' ~ 'M',
                               species == 'S_total_CFU' ~ 'S'),
           species_total_cfu = ifelse(is.numeric(species_total_cfu) & species_total_cfu == 0, NA, species_total_cfu)) %>% # if there are 0 CFU, replace with NA
    filter(rep != 6) %>% 
    ggplot(aes(x = factor(species_present, levels = c('S','SE', 'SM', 'SEM')), y = species_total_cfu))+
    stat_summary(fun = mean, shape = '-', size = 3, color = 'black')+
    geom_point(shape = 1)+
    theme_bw(12)+
    xlab("Treatment")+
    ylab("Total CFU")+
    scale_y_log10()+
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = 'none') +
    facet_grid(factor(species, levels = c('S','E', 'M'))~ factor(ecology, levels = c('Mutualism', 'Competition')))

saveRDS(object = total_pop_size_barseq, file = 'rds_plots/total_pop_size_barseq.rdata')


ggsave("plots/total_population_size_for_all_species.png",
       dpi = 300, width = 8, height = 4)



# there are significantly more cells in S cells in competitive S monouculture compared to competitive SM
df_total_long_target %>% 
    filter(rep != 6) %>% 
    mutate(
           species = case_when(species == 'E_total_CFU' ~ 'E',
                               species == 'M_total_CFU' ~ 'M',
                               species == 'S_total_CFU' ~ 'S')) %>% 
    filter(ecology == 'Competition') %>% 
    filter(species == 'S') %>% 
    filter(species_present %in% c('S','SM')) %>% 
    ungroup() %>% 
    t_test(formula = species_total_cfu ~ species_present)


df_total_long_target %>% 
    filter(rep != 6) %>% 
    mutate(
        species = case_when(species == 'E_total_CFU' ~ 'E',
                            species == 'M_total_CFU' ~ 'M',
                            species == 'S_total_CFU' ~ 'S')) %>% 
    filter(ecology == 'Competition') %>% 
    filter(species == 'S') %>% 
    filter(species_present %in% c('S','SM')) %>% 
    ungroup() %>% 
    ggplot(aes(species_present, species_total_cfu)) +
    geom_boxplot()
