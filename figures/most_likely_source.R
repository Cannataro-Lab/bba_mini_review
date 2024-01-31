



## find set of variants, highest prevalence, highest effect 

prev_cutoff <- 10
selection_top <- 10

luad_cesa$selection$selection.1 |> 
  arrange(desc(selection_intensity)) |> 
  slice_head(n = selection_top) -> 
  top_selected

cesa_chosen <-  luad_cesa$selection$selection.1 |> 
  filter(included_with_variant >= prev_cutoff | variant_name %in% top_selected$variant_name)


source_probs <- luad_attr$mutational_sources$source_probabilities |> 
  filter(variant_id %in% cesa_chosen$variant_id)


# most likely source from each tumor
source_probs |>  
  pivot_longer(cols = -c(variant_id,Unique_Patient_Identifier)) |> 
  filter(value>0) |> 
  ggplot(aes(x=variant_id, y=value,fill = name)) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle=90))



source_probs |>  
  pivot_longer(cols = -c(variant_id,Unique_Patient_Identifier)) |> 
  group_by(variant_id, Unique_Patient_Identifier) |> 
  slice_max(value) |> 
  ggplot(aes(x=variant_id,fill=name)) + 
  geom_bar()
