rm(list = ls())

library(tidyverse)

HPLC_path <- 'data/HPLC_data/'

meta <- read_csv(file = 'data/HPLC_data/spent_media_HPLC_submission_20230130.csv')

txt_files <- list.files(HPLC_path, pattern = '\\.txt$', full.names = TRUE)

# extract the first three characters of the file name to get get the tube ID
short_file_names <- str_sub(basename(txt_files), start = 1, end = 3)

# function that we will map to format all of the HPLC data
make_hplc_tibble <- function(file){
    #read the file as a character vector
    file_content <- read_lines(file)
    #remove the first 42 lines -- this contains meta information from the HPLC run
    file_content <- file_content[-(1:42)]
    #convert the character vector into a single string with newline separators
    file_string <- paste(file_content, collapse = '\n')
    #read the tab separated data into a tibble
    tibble_data <- read_tsv(file_string)
    return(tibble_data)
}

data_list <- map(txt_files, make_hplc_tibble)

data_list_named <- set_names(data_list, short_file_names)

# Combine all elements of data_list_named, adding a new column for the element names
df_all <- list_rbind(data_list_named, names_to = 'Tube_Number') %>% 
    left_join(meta, by = 'Tube_Number')

df_all %>% 
    filter(Tube_Number != 'WHC') %>% 
    filter(!Tube_Number %in% c('WH1', 'WH2', 'WH6')) %>% 
    # filter(Tube_Number %in% c('WH4', 'WH7')) %>% # filter out just the midlog
    ggplot(aes(x = `Time (min)`, y = `Value (mAU)`, color = Tube_Number))+
    geom_line() +
    # scale_y_log10() +
    # ylim(0, 5)+
    xlim(0, 5)+
    theme_bw(12)+ #+
    facet_wrap(~Description, ncol = 3)



df_all %>% 
    # filter(Tube_Number != 'WHC') %>% 
    # filter(!Tube_Number %in% c('WH1', 'WH2', 'WH6')) %>% 
    filter(Tube_Number %in% c('WH4', 'WH8')) %>% # filter out just the midlog
    ggplot(aes(x = `Time (min)`, y = `Value (mAU)`, color = Description))+
    geom_line() +
    # scale_y_log10() +
    # ylim(0, 5)+
    xlim(0, 5)+
    theme_bw(12) + 
    facet_wrap(~ Description, ncol = 1)


df_simple <- df_all %>% 
    filter(Tube_Number %in% c('WH4', 'WH8')) %>% 
    mutate(
        simple_descript = case_when(Tube_Number == 'WH4' ~ 'E spent media',
                                    Tube_Number == 'WH8' ~ 'S grown in E spent media'))

df_simple %>% 
    ggplot(aes(x = `Time (min)`, y = `Value (mAU)`))+
    geom_line() +
    xlim(0,5) +
    theme_bw(12) +
    facet_wrap(~ simple_descript, ncol = 1)

ggsave('plots/hplc_spent.png',
       dpi = 300, height = 2.6, width = 2.6)

