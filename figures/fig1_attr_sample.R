
library(patchwork)

dataset_name = "luad_cesa"
sample = "TCGA-73-4658"

source("output_data/plots/stacked_bar_plot.R")

fig1 <- tcga_to_bar_plot(dataset_name = "luad_cesa", sample = "TCGA-73-4658",
                         # highlight any contexts in attr plot
                         highlight_context = c("T[C>T]A","C[C>A]A"))

# 
# # outputs list of all SBS plots 
# fig1$sbs_plots

fig1$original_subs_plot + patchwork::wrap_plots(fig1$sbs_plots,ncol = 1) + fig1$attr_plot

ggsave(filename = "figures/fig1_scratch.pdf",width = 32,height = 8)


