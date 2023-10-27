# Mutualism reduces the severity of gene disruptions in predictable ways across microbial communities
![image](https://github.com/JonMartinson/ecology_DFE/assets/69863285/17420f3d-fc09-4e0b-87a6-2c3a62a5d833)





### Authors:
Jonathan N.V. Martinson☼, Jeremy M. Chacón☼, Brian A. Smith, Alex R. Villarreal, Ryan C. Hunter, William R. Harcombe

☼ contributed equally

---

## Repository Overview

This code repository contains the necessary data and scripts to reproduce the main figures and analyses presented in the manuscript, ["Mutualism Reduces the Severity of Gene Disruptions in Predictable Ways Across Microbial Communities."](https://www.nature.com/articles/s41396-023-01534-6)

Raw sequencing data is available in the NCBI BioProject: [PRJNA1008691](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1008691).

---

## Getting Started

1. Clone or download the repository to your local machine.
2. If you are using R studio, open the R project file, `DFE_analyses.Rproj`. This will set up your working directories correctly, so there's no need to change them manually (tested on macos/linux).
3. Unzip the `all.poolcount.zip` file in the `data` folder to access the data required for calculating gene fitnesses.

## BarSeq Fitness calculations
To genenerate fitnesses from the `all.poolcount` file (an output of the Feba pipeline described here doi.org/10.1128/mBio.00306-15), run the `mult_t0_pipeline.R` script located within the `scripts_2023/pipeline` folder. Fitness calculations are a hybrid of the methods used in doi.org/10.1128/mBio.00306-15 and doi.org/10.1038/s41564-020-00800-z. The resulting file `data/S0_barseq_fitnesses_20220722.csv` is ~ 2.26 Gb.

## Figures

The figures are saved as ggplot objects and pdfs in the 'rds_plots' and 'pdf_plots' folders, respectively. 

To recreate the figures, run the following scripts within the `scripts_2023` folder:

| Figure     | Script                                                      |
|------------|-------------------------------------------------------------|
| Figure 1A  | NA -- diagram                                             |
| Figure 1B,C| `01_barseq_dfe_analysis.R`                                  |
| Figure 2   | `02_barseq_scatterplots_barcharts.R`                        |
| Figure 3   | `03_bootstrap_model_proportion_stats.R`                     |
| Figure 4A  | `04_ORA_barseq.R`                                           |
| Figure 4B,C| `05_plot_individual_gene_fitness.R`                         |
| Figure 5A  | `05_plot_individual_gene_fitness.R`                         |
| Figure 5B  | `06_salmonella_mutant_ilvA_panC_spent_media_endpoint_OD.R`  |
| Figure 5C  | `07_ile_vitB5_gradient_analysis.R`                          |
| Figure 5D  | `08_Salmonella_spent_assay_ilvA_panC.R`                     |
| Figure S1  | NA -- diagram                                              |
| Figure S2  | `09_agar_scanner_growth.R` *                                 |
| Figure S3  | `10_wetlab_growth.R`                                        |
| Figure S4  | `11_median_fitness_effect.R`                                |
| Figure S5  | `12_heat_map_model_effect_on_genes.R`                       |
| Figure S6A | `13_S_aceA_plate_reader.R`                                  |
| Figure S6B | `14_hplc_plots.R`                                           |
| Figure S7  | `15_spent_media_pdxB_nadC_thiE.R`                           |
| Figure S8  | `16_mean_fitness_vs_mean_growth_rate.R`                     |

To save the plots as pdfs with consistent sizing use: `17_combined_plot_panels.R`

*Lawn density values calculated in `lawn_analysis_for_ecology_dfe_paper.ipynb`. The raw scanner images are available upon request. Additionally, we made a commandline python program for analyzing the scanner lawn densities: `scripts_2023/lawn_density_test`. `set_up_lawn_density.txt` contains instructions for use on a cluster or local terminal. `example_images/` contains a subset of the scanner images used in the experiment (note that the .tif files are compressed). 



---

## Required R Packages

Before running the scripts, ensure you have the following R packages installed:

- `tidyverse`: A collection of R packages designed for data science.
- `lubridate`: Functions to work with date-time data.
- `ggpubr`: Functions for creating publication-ready plots.
- `ggtext`: Enhances text rendering in `ggplot2` plots, allowing markdown formatting.
- `rstatix`: Functions for performing various statistical tests.
- `ggpmisc`: Miscellaneous functions for `ggplot2`.
- `clusterProfiler`: Tools for functional enrichment analysis.
- `ggridges`: Functions for creating ridgeline plots.
- `plater`: A package for working with plate-based data.


To install the required packages, run the following command in your R console:

```R
install.packages(c("tidyverse", "lubridate", "ggpubr", "rstatix", "ggpmisc", "clusterProfiler", "ggridges", "plater", "ggtext"))
