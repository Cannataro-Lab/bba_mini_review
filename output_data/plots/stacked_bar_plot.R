library(tidyverse)
library(plotly)


# dataset_name of the form luad_cesa, lusc_cesa etc.
tcga_to_bar_plot <- function(dataset_name, sample, highlight_context=NULL,subs_to_plot=c("C>A","C>G","C>T")){
  
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

  to_keep <- NULL
  
  for(sub_ind in 1:length(subs_to_plot)){
    to_keep <- c(to_keep, which(stringr::str_detect(pattern = subs_to_plot[sub_ind],
                                              string = trinuc_counts_samples_long$Type)))
    
  }
  
    
  trinuc_counts_samples_long_for_plot <- trinuc_counts_samples_long[unique(to_keep),]
    
  
  
  original_subs <- trinuc_counts_samples_long_for_plot |> 
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
  
  attr_sigs$Signatures <- factor(attr_sigs$Signatures, levels = signatures$Signatures)
  
  to_keep <- NULL
  
  for(sub_ind in 1:length(subs_to_plot)){
    to_keep <- c(to_keep, which(stringr::str_detect(pattern = subs_to_plot[sub_ind],
                                                    string = attr_sigs$Type)))
    
  }
  
  attr_sigs_for_plot <- attr_sigs[unique(to_keep),]
  
  
  if(!is.null(highlight_context)){
    
    # highlight_context <- c("T[C>T]A")
    # attr_sigs$bordered <- NA
    attr_sigs_for_plot <- attr_sigs_for_plot |> 
      mutate(bordered = case_when(as.character(Type) %in% highlight_context ~ "highlight", 
             TRUE ~ "no"))
    
    
    attr_plot <- ggplot(attr_sigs_for_plot)   
    attr_plot <- attr_plot + 
      geom_bar(aes(x=Type, 
                   y=absolute_signature, 
                   fill = Signatures,color=bordered), stat = "identity", position = "stack") + 
      scale_color_manual(values = c("black","white")) + 
      guides(colour = "none")
    
  }else{
    
    attr_plot <- ggplot(attr_sigs_for_plot) 
    attr_plot <- attr_plot + 
      geom_bar(aes(x=Type, y=absolute_signature, fill = Signatures), stat = "identity", position = "stack") 
    
  }
  
  attr_plot <- attr_plot + 
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
  
  to_keep <- NULL
  
  for(sub_ind in 1:length(subs_to_plot)){
    to_keep <- c(to_keep, which(stringr::str_detect(pattern = subs_to_plot[sub_ind],
                                                    string = this_sig_dist$Type)))
    
  }
  
  this_sig_dist <- this_sig_dist[unique(to_keep),]
  
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
