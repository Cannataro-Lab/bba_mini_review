
library(patchwork)

# dataset_name = "luad_cesa"
# sample = "TCGA-73-4658"

source("output_data/plots/stacked_bar_plot.R")

# fig1 <- tcga_to_bar_plot(dataset_name = "luad_cesa", sample = "TCGA-73-4658",
#                          # highlight any contexts in attr plot
#                          highlight_context = c("T[C>T]A","C[C>A]A"),subs_to_plot = c("C>A","C>T"))
# 
# # 
# # # outputs list of all SBS plots 
# # fig1$sbs_plots
# 
# fig1$original_subs_plot + patchwork::wrap_plots(fig1$sbs_plots,ncol = 1) + fig1$attr_plot
# 
# ggsave(filename = "figures/fig1_scratch.pdf",width = 20,height = 8)
# 
# 

# find LUAD w/ KRAS G12C for explanatory figure 

luad_kras_g12c <- luad_cesa$maf |> 
  filter(top_consequence == "KRAS_G12C")  |> 
  pull(Unique_Patient_Identifier)


KRAS_snvs <- luad_cesa@mutations$amino_acid_change |> 
  filter(gene == "KRAS") |> 
  filter(aachange == "G12C")

kras_g12c_trinuc_context <- luad_cesa@mutations$snv |> 
  filter(snv_id == KRAS_snvs$constituent_snvs[[1]]) |> 
  pull(trinuc_mut)

for(plot_ind in 1:length(luad_kras_g12c)){
  
  this_sample <- luad_kras_g12c[plot_ind]
  
  
  
  fig1 <- tcga_to_bar_plot(dataset_name = "luad_cesa", 
                           sample = this_sample,
                           # highlight any contexts in attr plot
                           highlight_context = c(kras_g12c_trinuc_context),
                           subs_to_plot = c("C>A","C>T"))
  
  # 
  # # outputs list of all SBS plots 
  # fig1$sbs_plots
  
  fig1$original_subs_plot + patchwork::wrap_plots(fig1$sbs_plots,ncol = 1) + fig1$attr_plot
  
  # save w specific name 
  ggsave(filename = paste0("figures/fig_1_tests/",this_sample,".pdf"),
         width = 20,height = 8)
  
  
}

# 75-5126 looks good for a figure without being too busy


source("output_data/plots/stacked_bar_plot.R")

fig1 <- tcga_to_bar_plot(dataset_name = "luad_cesa", 
                         sample = "TCGA-75-5126",
                         # highlight any contexts in attr plot
                         highlight_context = c(kras_g12c_trinuc_context),
                         subs_to_plot = c("C>A","C>T"))

fig1b <- patchwork::wrap_plots(fig1$sbs_plots$SBS2 + labs(tag = "B"),
                               fig1$sbs_plots$SBS4,
                               fig1$sbs_plots$SBS5,
                               fig1$sbs_plots$SBS13,ncol = 1)

fig1$original_subs_plot + labs(tag = "A") + 
  fig1b + 
  fig1$attr_plot + labs(tag = "C")

ggsave(filename = "figures/fig1.pdf",width = 15,height = 8)
ggsave(filename = "figures/fig1.png",width = 15,height = 8)

luad_cesa$maf |> 
  filter(Unique_Patient_Identifier == "TCGA-75-5126") |> 
  filter(top_gene == "KRAS")


luad_cesa$maf |> 
  filter(Unique_Patient_Identifier == "TCGA-75-5126") |> 
  left_join(luad_cesa$selection$selection.1, by = c("top_consequence" = "variant_name")) |> 
  filter(!is.na(selection_intensity)) |> 
  arrange(desc(selection_intensity))


luad_cesa$maf |> 
  filter(Unique_Patient_Identifier == "TCGA-75-5126") |> 
  nrow()




luad_cesa$selection$selection.1 |>  
  filter(variant_name == "MTOR_R2266P")
