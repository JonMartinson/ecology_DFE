# Mutualism Reduces the Severity of Gene Disruptions in Predictable Ways Across Microbial Communities
![image](https://user-images.githubusercontent.com/69863285/235250401-e4a04097-b5f6-4edf-bba0-e56a6d054f4d.png)
### Authors:
Jonathan N.V. Martinson, Jeremy M. Chacón, Brian A. Smith, Alex R. Villarreal, Ryan C. Hunter, William R. Harcombe

---

## Repository Overview

This code repository contains the necessary data and scripts to reproduce the main figures and analyses presented in the manuscript, "Mutualism Reduces the Severity of Gene Disruptions in Predictable Ways Across Microbial Communities." To reproduce the analyses, clone or download the project and follow the instructions provided below.

---

## Getting Started

1. Clone or download the repository to your local machine.
2. If you are using R studio, open the R project file, `DFE_analyses.Rproj`. This will set up your working directories correctly, so there's no need to change them manually.
3. Unzip the `all.poolcount.zip` file in the `data` folder to access the data required for calculating gene fitnesses.

## Fitness calculation
To genenerate fitnesses from the all.poolcount file (an output of the Feba pipeline described here http://mbio.asm.org/content/6/3/e00306-15.full), run the 'mult_t0_pipeline.R' script located within the 'scripts_2023/pipeline' folder.

## Figures

The figures are saved as gg objects and pdfs in the 'rds_plots' and 'pdf_plots' folders, respectively. 

## Recreating Figures

To recreate the figures, run the following scripts within the `scripts_2023` folder:

| Figure      | Script                                             |
|-------------|----------------------------------------------------|
| Figure 1B,C | `barseq_dfe_analysis.R`                            |
| Figure 2A   | `barseq_scatterplots_label_change.R`               |
| Figure 2B   | `barseq_gene_effects_upset.R`                      |
| Figure 3A   | `bootstrap_model_proportion_stats.R`               |
| Figure 3B   | `barseq_dfe_analysis.R`                            |
| Figure 4    | `ORA_barseq.R`                                     |
| Figure 5    | `plot_individual_gene_fitness.R`                   |
| Figure 6A   | `plot_individual_gene_fitness.R`                   |
| Figure 6B   | `S_aceA_plate_reader_figure_modified.R`            |
| Figure 6C   | `hplc_plots.R`                                     |
| Figure 7A   | `plot_individual_gene_fitness.R`                   |
| Figure 7B   | `salmonella_mutant_ilvA_panC_spent_media_endpoint_OD.R` |
| Figure 7C   | `ile_vitB5_gradient_analysis.R`                    |
| Figure 7D   | `spent_assay_ilvA_panC_20230417.R`                 |


---

## Required R Packages

Before running the scripts, ensure you have the following R packages installed:

- `tidyverse`: A collection of R packages designed for data science.
- `lubridate`: Functions to work with date-time data.
- `ggpubr`: Functions for creating publication-ready plots.
- `rstatix`: Functions for performing various statistical tests.
- `ggpmisc`: Miscellaneous functions for `ggplot2`.
- `clusterProfiler`: Tools for functional enrichment analysis.
- `ggridges`: Functions for creating ridgeline plots.
- `plater`: A package for working with plate-based data.

To install the required packages, run the following command in your R console:

```R
install.packages(c("tidyverse", "lubridate", "ggpubr", "rstatix", "ggpmisc", "clusterProfiler", "ggridges", "plater"))
