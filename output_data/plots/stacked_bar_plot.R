library(tidyverse)
library(plotly)

luad_attr$effect_shares

first_sample <- luad_cesa$trinuc_rates %>%
  filter(Unique_Patient_Identifier == "TCGA-05-4244")

signatures <- read.table("COSMIC_v3.4_SBS_GRCh37.txt", header = TRUE)
signatures_longer <- signatures %>%
  pivot_longer(cols = starts_with("SBS"), names_to = "Signatures", values_to = "raw_signature")

signatures_longer %>%
  filter(Signatures == "SBS1") %>%
  filter(Values > 0.1) %>%
  ggplot() +
  geom_col(aes(x=Type, y=Values))


first_example <- luad_cesa@trinucleotide_mutation_weights[["raw_signature_weights"]] %>%
  filter(Unique_Patient_Identifier == "TCGA-05-4244") %>%
  pivot_longer(cols = starts_with("SBS"), names_to = "Signatures", values_to = "raw_signature_weights")

highest_sigs <- signatures_longer %>%
  right_join(first_example, by = "Signatures") %>%
  mutate(absolute_signature = raw_signature*raw_signature_weights) %>%
  group_by(Type) %>%
  summarize(sum = sum(absolute_signature)) %>%
  filter(sum > 4)

plot_try <- signatures_longer %>%
  right_join(first_example, by = "Signatures") %>%
  mutate(absolute_signature = raw_signature*raw_signature_weights) %>%
  filter(Type %in% highest_sigs$Type) %>%
  ggplot() +
  geom_bar(aes(x=Type, y=absolute_signature, fill = Signatures), stat = "identity", position = "stack")

ggplotly(plot_try)
