![image](https://user-images.githubusercontent.com/69863285/235230862-9966ba3a-8389-4459-9ddf-47440e99c1b7.png)




# Mutualism Reduces the Severity of Gene Disruptions in Predictable Ways Across Microbial Communities

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
