library(tidyverse)
library(plotly)


# dataset_name of the form luad_cesa, lusc_cesa etc.
tcga_to_bar_plot <- function(dataset_name, sample){
  
  dataset <- get(dataset_name)
  
  signatures <- read.table("COSMIC_v3.4_SBS_GRCh37.txt", header = TRUE)
  signatures_longer <- signatures %>%
    pivot_longer(cols = starts_with("SBS"), names_to = "Signatures", values_to = "raw_signature")
  
  sample_example <- dataset@trinucleotide_mutation_weights[["raw_signature_weights"]] %>%
    filter(Unique_Patient_Identifier == sample) %>%
    pivot_longer(cols = starts_with("SBS"), names_to = "Signatures", values_to = "raw_signature_weights")
  
  sample_highest_sigs <- signatures_longer %>%
    right_join(sample_example, by = "Signatures") %>%
    mutate(absolute_signature = raw_signature*raw_signature_weights) %>%
    group_by(Type) %>%
    summarize(sum = sum(absolute_signature)) %>%
    slice_max(sum, n = 10)
  
  signatures_longer %>%
    right_join(sample_example, by = "Signatures") %>%
    mutate(absolute_signature = raw_signature*raw_signature_weights) %>%
    filter(Type %in% sample_highest_sigs$Type) %>%
    ggplot() +
    geom_bar(aes(x=Type, y=absolute_signature, fill = Signatures), stat = "identity", position = "stack") +
    labs(x="Number of substitutions", y="Trinucleotide context")
}
