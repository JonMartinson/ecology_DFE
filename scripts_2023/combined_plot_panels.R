library(tidyverse)
library(ggpubr)

# Figure 1

# fig 1A was made in powerpoint
fig_1B <- readRDS(file = 'rds_plots/DFE_histogram.rdata') 
fig_1C <- readRDS(file = 'rds_plots/mean_fit_effect_plot.rdata')

ggsave(fig_1B, filename = 'pdf_plots/fig_1b.pdf', device = 'pdf', width = 6, height = 3)
ggsave(fig_1C, filename = 'pdf_plots/fig_1c.pdf', device = 'pdf', width = 6, height = 3)

 # Figure 2


fig_2ABDE <- readRDS(file = 'rds_plots/scatterplot_cocult_vs_mono.rdata')
fig_2CF<- readRDS( file = 'rds_plots/sig_gene_overlap.rdata')

ggsave(fig_2ABDE, filename = 'pdf_plots/fig_2_ABDE.pdf', device = 'pdf', width = 6, height = 6)
ggsave(fig_2CF, filename = 'pdf_plots/fig_2CF.pdf', device = 'pdf', width = 3, height = 6)

# Figure 3

# fig_3A_i <-  readRDS(file = 'rds_plots/mut_box_bootstrap.rdata')
# fig_3A_ii <-  readRDS(file = 'rds_plots/comp_box_bootstrap.rdata') 
fig_3A <- readRDS(file = 'rds_plots/bootstrap_models_combined.rdata')
fig_3B <- readRDS(file = 'rds_plots/additivity_bar_plot.rdata')

ggsave(fig_3A, filename = 'pdf_plots/fig_3A.pdf', device = 'pdf', height = 2, width = 5)
ggsave(fig_3B, filename = 'pdf_plots/fig_3B.pdf', device = 'pdf', height = 2, width = 4)
# Figure 4

fig_4 <- readRDS(file = 'rds_plots/mutualism_ORA.rdata')

ggsave(fig_4, filename = 'pdf_plots/fig_4.pdf', device = 'pdf', height = 4, width = 8)
# Figure 5

fig_5 <- readRDS(file = 'rds_plots/nitrogen_barseq.rdata')

ggsave(fig_5, filename = 'pdf_plots/fig_5.pdf', device = 'pdf', height = 4, width = 4)

# Figure 6


fig_6A <- readRDS(file = 'rds_plots/ace_gal_barseq.rdata')
fig_6B <-  readRDS('rds_plots/aceA_growth_curves.rdata')
fig_6C<- readRDS(file = 'rds_plots/hplc_simple.rdata')+ labs(x = 'Retention Time (min)', y = 'Value (mAU v210nm)')

fig_6_all<- ggarrange(fig_6A,
          fig_6B,
          fig_6C, 
          ncol = 3,
          align = 'h',widths = c(1.5, 2, 1.25))

ggsave(fig_6_all, filename = 'pdf_plots/fig_6.pdf', device = 'pdf', height = 4, width = 10)

# Figure 7 
fig_7A <- read_rds(file = 'rds_plots/simple_amino_vit_barseq.rdata') + labs(y = 'Fitness Normalized') 
fig_7B <- read_rds(file = 'rds_plots/EM_spent_ilvA_panC.rdata') 
fig_7C_i <- read_rds(file = 'rds_plots/isoleucine_gradient.rdata')
fig_7C_ii <- read_rds(file = 'rds_plots/vitamin_B5_gradient.rdata') 
fig_7D <- read_rds(file = 'rds_plots/S_spent_media_ilvA_panC.rdata')


fig_7C <- ggarrange(fig_7C_i,fig_7C_ii + rremove("ylab") +theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank()),common.legend = TRUE,
    widths = c(1, .87))

fig_7B_to_7D<- ggarrange(fig_7B, fig_7C, fig_7D, ncol = 1)

fig_7 <- ggarrange(fig_7A, fig_7B_to_7D, widths = c(.8, 1.75))

ggsave(fig_7, filename = 'pdf_plots/fig_7.png', dpi = 1200, width = 9, height = 7, device = 'png')
# DID NOT SAVE AS PDF BECAUSE GREEK SYMBOL ∆ DOES NOT RENDER.

# ggsave(fig_7, filename = 'pdf_plots/fig_7.png', dpi = 1200, width = 9, height = 9, device = 'png')
# ggsave(fig_7, filename = 'pdf_plots/fig_7.pdf', width = 9, height = 9, device = 'pdf')
# ggsave( filename = 'pdf_plots/fig_7.pdf', dpi = 1200, device = cairo_pdf)


 # Supp Fig 1

Sfig_1A <-  readRDS(file = 'rds_plots/total_pop_size_barseq.rdata') 
Sfig_1B <- readRDS(file = 'rds_plots/generations_species.rdata')
Sfig_1C <- readRDS(file= 'rds_plots/species_frequencies_end.rdata')
Sfig_1AB<- ggarrange(Sfig_1A,Sfig_1B, ncol = 2)
Sfig_1 <- ggarrange(Sfig_1AB,Sfig_1C, ncol = 1)

ggsave(Sfig_1, filename = 'pdf_plots/Sfig_1.pdf', device = 'pdf', height = 8, width = 7)

# Supp Fig 2

Sfig_2A <- readRDS(file = 'rds_plots/lawn_growthcurve.rdata')
Sfig_2B <- readRDS(file = 'rds_plots/lawn_growth_rate.rdata')

Sfig_2 <- ggarrange(Sfig_2A, Sfig_2B)

ggsave(plot = Sfig_2, file = 'pdf_plots/Sfig_2.pdf', height = 3, width = 8, device = 'pdf' )

# Supp Fig 3 
Sfig_3A<- readRDS('rds_plots/mutualism_gene_heatmap.rdata')
Sfig_3B<- readRDS('rds_plots/competition_gene_heatmap.rdata')

ggsave(plot = Sfig_3A, file = 'pdf_plots/Sfig_3A.pdf', height = 12, width = 3, device = 'pdf' )
ggsave(plot = Sfig_3B, file = 'pdf_plots/Sfig_3B.pdf', height = 12, width = 3, device = 'pdf' )

# Supp Fig 4 

Sfig_4 <- readRDS('rds_plots/ilvA_panC_supplemented_yield.rdata')


ggsave(plot = Sfig_4, file = 'pdf_plots/Sfig_4.png', height = 2, width = 4, device = 'png', dpi = 1200 )
# ggsave(plot = Sfig_4, file = 'pdf_plots/Sfig_4.pdf', height = 3, width = 5, device = 'pdf' )
