library(tidyverse)
library(ggpubr) # arranging figures in multipanel orientation
library(ggtext) # allows markdown text formatting in figures

#### Main Text Figures ####

# Figure 1

# fig 1A was made in powerpoint
fig_1B <- readRDS(file = 'rds_plots/DFE_histogram.rdata')
fig_1C <- readRDS(file = 'rds_plots/mean_fit_effect_plot.rdata')

ggsave(fig_1B, filename = 'pdf_plots/fig_1b.pdf', device = 'pdf', width = 6, height = 3)
ggsave(fig_1C, filename = 'pdf_plots/fig_1c.pdf', device = 'pdf', width = 6, height = 3)

 # Figure 2

fig_2ABDE <- readRDS(file = 'rds_plots/scatterplot_cocult_vs_mono.rdata')
fig_2CF<- readRDS( file = 'rds_plots/fig_2_bar.rdata')

ggsave(fig_2ABDE, filename = 'pdf_plots/fig_2_ABDE.pdf', device = 'pdf', width = 6, height = 6)
ggsave(fig_2CF, filename = 'pdf_plots/fig_2CF.pdf', device = 'pdf', width = 5, height = 6)

# Figure 3

fig_3A <- readRDS(file = 'rds_plots/bootstrap_models_combined.rdata')
fig_3B <- readRDS(file = 'rds_plots/additivity_bar_plot.rdata')

ggsave(fig_3A, filename = 'pdf_plots/fig_3A.pdf', device = 'pdf', height = 2, width = 5)
ggsave(fig_3B, filename = 'pdf_plots/fig_3B.pdf', device = 'pdf', height = 2, width = 4)
# Figure 4

fig_4a <- readRDS(file = 'rds_plots/mutualism_ORA.rdata')
fig_4BC <- readRDS('rds_plots/nitro_galace_barseq_gene_fit.rdata')

ggsave(fig_4a, filename = 'pdf_plots/fig_4a.pdf', device = 'pdf', height = 7, width = 8)
ggsave(fig_4BC, filename = 'pdf_plots/fig_4BC.pdf', device = 'pdf', height = 7, width = 4)


# Figure 5


fig_5A <- read_rds(file = 'rds_plots/simple_amino_vit_barseq.rdata') + labs(y = 'Fitness effect') 
fig_5B <- read_rds(file = 'rds_plots/EM_spent_ilvA_panC.rdata') 
fig_5C_i <- read_rds(file = 'rds_plots/isoleucine_gradient.rdata')
fig_5C_ii <- read_rds(file = 'rds_plots/vitamin_B5_gradient.rdata') 
fig_5D <- read_rds(file = 'rds_plots/S_spent_media_ilvA_panC.rdata')


fig_5C <- ggarrange(fig_5C_i,fig_5C_ii + rremove("ylab") +theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank()),common.legend = TRUE,
    widths = c(1, .87))

fig_5B_to_5D<- ggarrange(fig_5B, fig_5C, fig_5D, ncol = 1)

fig_5 <- ggarrange(fig_5A, fig_5B_to_5D, widths = c(1, 1.75))

ggsave(fig_5, filename = 'pdf_plots/fig_5.png', dpi = 1200, width = 10, height = 7, device = 'png')

#### Supplementary Figures ####

# Supplementary Figure S2

fig_S2A <- readRDS(file = 'rds_plots/lawn_growthcurve.rdata')
fig_S2B <- readRDS(file = 'rds_plots/lawn_growth_rate.rdata')

fig_S2 <- ggarrange(fig_S2A, fig_S2B)

ggsave(plot = fig_S2, file = 'pdf_plots/fig_S2.pdf', height = 3, width = 8, device = 'pdf' )


# Supplementary Figure S3

fig_S3A <-  readRDS(file = 'rds_plots/total_pop_size_barseq.rdata') 
fig_S3B <- readRDS(file = 'rds_plots/generations_species.rdata')
fig_S3C <- readRDS(file= 'rds_plots/species_frequencies_end.rdata')
fig_S3AB<- ggarrange(fig_S3A,fig_S3B, ncol = 2)
fig_S3 <- ggarrange(fig_S3AB,fig_S3C, ncol = 1)

ggsave(fig_S3, filename = 'pdf_plots/fig_S3.pdf', device = 'pdf', height = 8, width = 9)

    

# Supplementary Figure S4
# medians close to 0

fig_S4 <-  readRDS(file = 'rds_plots/median_fitness_effect.rdata') 
ggsave(fig_S4, filename = 'pdf_plots/fig_S4.pdf', device = 'pdf', height = 3, width = 5)

# Supplementary Figure S5

fig_S5A <- readRDS('rds_plots/mutualism_gene_heatmap.rdata')
fig_S5B <- readRDS('rds_plots/competition_gene_heatmap.rdata')

ggsave(plot = fig_S5A, file = 'pdf_plots/fig_S5A.pdf', height = 12, width = 16, device = 'pdf' )
ggsave(plot = fig_S5B, file = 'pdf_plots/fig_S5B.pdf', height = 12, width = 9, device = 'pdf' )

# Supplementary Figure S6
fig_S6A <- readRDS('rds_plots/aceA_growth_curves.rdata')
fig_S6B <- readRDS('rds_plots/hplc_simple.rdata')

fig_S6 <- ggarrange(fig_S6A,fig_S6B, ncol = 1) 

ggsave(plot = fig_S6, file = 'pdf_plots/fig_S6.png', height = 6.5, width = 9, device = 'png', units = 'in', dpi = 1200)

# Supplementary Figure S7

fig_S7 <- readRDS('rds_plots/BEI_spent_media.rdata')

ggsave(plot = fig_S7, file = 'pdf_plots/fig_S7.png', height = 4, width = 9, device = 'png', units = 'in', dpi = 1200)

# Supplementary Figure S8

fig_S8 <- readRDS('rds_plots/growth_vs_fitness.rdata')

ggsave(plot = fig_S8, file = 'pdf_plots/fig_S8.pdf', height = 2.5, width = 5, device = 'pdf')
