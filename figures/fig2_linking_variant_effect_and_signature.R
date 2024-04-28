


# generate figure 2 linking high effect/prevalent variants with 
# mutational signatures and likely mutagenic sources




luad_attr_plot <- sources_of_variants(cesa_dataset_name = "luad_cesa",
                                      cesa_attr = luad_attr,selection_top = 12,
                                      prev_cutoff = 10) + 
  labs(title = "Lung adenocarcinoma")

lihc_attr_plot <- sources_of_variants(cesa_dataset_name = "lihc_cesa",
                                      cesa_attr = lihc_attr,selection_top = 11,
                                      prev_cutoff = 10) + 
  labs(title = "Liver Hepatocellular Carcinoma")

luad_attr_plot / lihc_attr_plot + plot_annotation(tag_levels = "A")

ggsave(filename = "figures/fig2.png",width = 10,height = 6)
ggsave(filename = "figures/fig2.tif",width = 10,height = 6)




