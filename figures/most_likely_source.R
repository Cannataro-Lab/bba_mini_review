



## find set of variants, highest prevalence, highest effect 


# cesa_dataset_name <- "luad_cesa"
# cesa_attr <- luad_attr



sources_of_variants <- function(cesa_dataset_name, cesa_attr, 
                                prev_cutoff = 10, 
                                selection_top = 10){
  
  
  dataset <- get(cesa_dataset_name)
  
  ordered_data <- dataset$selection$selection.1 |>
    arrange(desc(selection_intensity))
  
  ordered_data$variant_plot_name <- paste(stringr::str_split_i(string = ordered_data$variant_id,pattern = "_",i = 1),
                                               stringr::str_split_i(string = ordered_data$variant_id,pattern = "_",i = 2))
  
  
  dataset$selection$selection.1 |>
    arrange(desc(selection_intensity)) |> 
    slice_head(n = selection_top) -> 
    top_selected
  
  cesa_chosen <-  dataset$selection$selection.1 |> 
    filter(included_with_variant >= prev_cutoff | variant_name %in% top_selected$variant_name)
  
  
  source_probs <- cesa_attr$mutational_sources$source_probabilities |> 
    filter(variant_id %in% cesa_chosen$variant_id)
  
  
  source_probs_long <- source_probs |> 
    pivot_longer(cols = -c(variant_id,Unique_Patient_Identifier)) |> 
    filter(value>0)
  
  
  source_probs_long
  
  source_probs_long <- source_probs_long |> 
    mutate(signature = name) |> 
    mutate(signature_process = case_when(
      signature == "SBS1" ~ "Deamination with age, clock-like (1)",
      signature == "SBS5" ~ "Unknown, clock-like (5)",
      signature %in% c("SBS2","SBS13") ~ "APOBEC (2,13)",
      signature %in% c("SBS3") ~ "Defective homologous recombination (3)",
      signature %in% c("SBS4","SBS29") ~ "Tobacco (4,29)",
      signature %in% c("SBS7a","SBS7b","SBS7c","SBS7d","SBS38") ~ "UV light (7a–d,38)",
      signature %in% c("SBS11","SBS31","SBS32","SBS35") ~ "Prior treatment (11,31,32,35)",
      signature %in% c("SBS22","SBS24","SBS42","SBS88") ~ "Mutagenic chemical exposure (22,24,42,88)",
      signature %in% c("SBS16") ~ "Alcohol-associated (16)",
      signature %in% c("SBS18") ~ "Damage by reactive oxygen species (18)", 
      TRUE ~ "Other signatures"))
  
  
  source_probs_long
  
  source_probs_long$variant_plot_name <- paste(stringr::str_split_i(string = source_probs_long$variant_id,pattern = "_",i = 1),
                                               stringr::str_split_i(string = source_probs_long$variant_id,pattern = "_",i = 2))
  
  source_probs_long$variant_plot_name <- factor(source_probs_long$variant_plot_name,
                                                levels = ordered_data$variant_plot_name)
  
  source_probs_long_summary <- source_probs_long |> 
    group_by(name, signature_process,variant_plot_name) |>  
    summarize(total_val = sum(value))
  
  
  source_probs_long_summary$signature_process <- 
    factor(source_probs_long_summary$signature_process,
           levels = names(color_vec))
  
  # most likely source from each tumor
  source_probs_long_summary |> 
    ggplot(aes(x=variant_plot_name, y=total_val,fill = signature_process),size=0) + 
    geom_col() + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle=90)) + 
    scale_fill_manual(values = color_vec) + 
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) + 
    labs(y="Number of substitutions",x="Variant", fill = "Signature process")
  
  
  
  
}

