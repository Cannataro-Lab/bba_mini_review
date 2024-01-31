library(tidyverse)
library(plotly)


# dataset_name of the form luad_cesa, lusc_cesa etc.
tcga_to_bar_plot <- function(dataset_name, sample){
  
  dataset <- get(dataset_name)
  
  trinuc_order <- colnames(ces.refset.hg38::ces.refset.hg38$signatures$COSMIC_v3.2$signatures)
  
  # signatures <- read.table("COSMIC_v3.4_SBS_GRCh37.txt", header = TRUE)
  signatures <- ces.refset.hg38::ces.refset.hg38$signatures$COSMIC_v3.2$signatures
  signatures <- signatures |> mutate(Signatures = rownames(signatures))
  
  signatures_longer <- signatures %>%
    pivot_longer(cols = -"Signatures", values_to = "raw_signature",names_to = "Type")
  
  
  trinuc_counts_samples <- dataset@trinucleotide_mutation_weights$trinuc_snv_counts
  
  count_types <- rownames(trinuc_counts_samples )
  
  trinuc_counts_samples <- as.data.frame(trinuc_counts_samples)
  trinuc_counts_samples$Type <- count_types
  
  trinuc_counts_samples_long <- trinuc_counts_samples |> 
    pivot_longer(cols = -Type, names_to = "Sample",values_to = "Number of substitutions")
  
  trinuc_counts_samples_long$Type <- factor(trinuc_counts_samples_long$Type, levels = trinuc_order) 
  
  original_subs <- trinuc_counts_samples_long |> 
    filter(Sample == sample) |> 
    ggplot(aes(x=Type, y=`Number of substitutions`)) + 
    geom_col() + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90)) + 
    labs(x="Trinucleotide context")
  
  
  
  sample_sbs_prop <- dataset$mutational_signatures$biological_weights |> 
    filter(Unique_Patient_Identifier == sample) |> 
    select(Unique_Patient_Identifier, starts_with("SBS")) |> 
    pivot_longer(-Unique_Patient_Identifier, names_to = "Signature",values_to = "Proportion") |> 
    filter(Proportion > 0)
  
  
  # sample_sbs_prop
  
  
  sample_example <- dataset$mutational_signatures$raw_attributions %>%
    filter(Unique_Patient_Identifier == sample) %>%
    pivot_longer(cols = starts_with("SBS"), names_to = "Signatures", values_to = "raw_signature_weights")
  
  sample_highest_sigs <- signatures_longer %>%
    right_join(sample_example, by = "Signatures") %>%
    mutate(absolute_signature = raw_signature*raw_signature_weights) %>%
    group_by(Type) %>%
    summarize(sum = sum(absolute_signature))
  
  # %>%
  #   slice_max(sum, n = 10)
  # 
  attr_sigs <- signatures_longer %>%
    right_join(sample_example, by = "Signatures") %>%
    mutate(absolute_signature = raw_signature*raw_signature_weights) %>%
    filter(Type %in% sample_highest_sigs$Type) %>%
    filter(absolute_signature >0) 
  
  attr_sigs$Type <- factor(attr_sigs$Type, levels = trinuc_order) 
  
  
  attr_plot <- ggplot(attr_sigs) +
    geom_bar(aes(x=Type, y=absolute_signature, fill = Signatures), stat = "identity", position = "stack") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90)) + 
    labs(y="Number of substitutions", x="Trinucleotide context") + 
    theme(legend.position=c(.9,.75)) + 
    scale_fill_manual(values = color_vec_sbs)
  
  
  
  # percent each sig plots
  
  sbs_plots <- vector(mode = "list",length = nrow(sample_sbs_prop))

  names(sbs_plots) <- sample_sbs_prop$Signature
  
  for(this_sbs_ind in 1:length(sbs_plots)){
    
    
  this_sbs <- sample_sbs_prop[this_sbs_ind,"Signature"] |> pull()
  this_prop <- sample_sbs_prop[this_sbs_ind,"Proportion"] |> pull()
  
  this_sig_dist <- signatures_longer |> 
    filter(Signatures == this_sbs)
  
  
  this_sig_dist$Type <- factor(this_sig_dist$Type, levels = trinuc_order) 
  
  
 
   sbs_ggplot <- ggplot(this_sig_dist) +
    geom_col(aes(x=Type, y=raw_signature), fill = color_vec_sbs[this_sbs]) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90)) + 
    labs(y="Number of substitutions", x="Trinucleotide context",
         title = this_sbs,
         subtitle = paste("Proportion signature attributed:", round(this_prop,3)))+ 
    theme(legend.position="none") 
   
   sbs_plots[[this_sbs_ind]] <- sbs_ggplot
  
  }
  
  
  return(
    list(
      original_subs_plot = original_subs,
      attr_plot = attr_plot,
      sbs_plots = sbs_plots
    )
  )
  
}
