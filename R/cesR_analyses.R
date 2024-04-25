# run cesR on MAF files from TCGA 
# (previously ran R/download_data.R)


library(cancereffectsizeR)
library(ces.refset.hg38)


# THCA ---- 


thca_maf <- preload_maf(maf = "output_data/THCA.maf",refset = "ces.refset.hg38")
  
cesa <- CESAnalysis(refset = "ces.refset.hg38") |> 
  load_maf(maf = thca_maf) |> 
  trinuc_mutation_rates(signature_set = "COSMIC_v3.2",
                        signature_exclusions = 
                          suggest_cosmic_signature_exclusions(cancer_type = "THCA", treatment_naive = TRUE),
                        cores = 6) |> 
  gene_mutation_rates(covariates = ces.refset.hg38$covariates$THCA) 


  # changed to be all, not just recurrent 
cesa <- cesa |> ces_variant(cores = 6)

cesa$selection$selection.1 |>  
  arrange(desc(selection_intensity))


save_cesa(cesa = cesa, file = "output_data/cesa_files/THCA_cesa.rds")


# LUAD ---- 

luad_maf <- preload_maf(maf = "output_data/LUAD.maf",refset = "ces.refset.hg38")

cesa <- CESAnalysis(refset = "ces.refset.hg38") |> 
  load_maf(maf = luad_maf) |> 
  trinuc_mutation_rates(signature_set = "COSMIC_v3.2",
                        signature_exclusions = 
                          suggest_cosmic_signature_exclusions(cancer_type = "LUAD", treatment_naive = TRUE),
                        cores = 6) |> 
  gene_mutation_rates(covariates = ces.refset.hg38$covariates$lung) |> 
  ces_variant(cores = 6)

cesa$selection$selection.1 |>  
  arrange(desc(selection_intensity))

save_cesa(cesa = cesa, file = "output_data/cesa_files/LUAD_cesa.rds")




# LUSC ---- 

lusc_maf <- preload_maf(maf = "output_data/LUSC.maf",refset = "ces.refset.hg38")

cesa <- CESAnalysis(refset = "ces.refset.hg38") |> 
  load_maf(maf = lusc_maf) |> 
  trinuc_mutation_rates(signature_set = "COSMIC_v3.2",
                        signature_exclusions = 
                          suggest_cosmic_signature_exclusions(cancer_type = "LUSC", treatment_naive = TRUE),
                        cores = 6) |> 
  gene_mutation_rates(covariates = ces.refset.hg38$covariates$lung) |> 
  ces_variant(cores = 6)

cesa$selection$selection.1 |>  
  arrange(desc(selection_intensity))

save_cesa(cesa = cesa, file = "output_data/cesa_files/LUSC_cesa.rds")



# ESCA ---- 

esca_maf <- preload_maf(maf = "output_data/ESCA.maf",refset = "ces.refset.hg38")

cesa <- CESAnalysis(refset = "ces.refset.hg38") |> 
  load_maf(maf = esca_maf) |> 
  trinuc_mutation_rates(signature_set = "COSMIC_v3.2",
                        signature_exclusions = 
                          suggest_cosmic_signature_exclusions(cancer_type = "Eso-AdenoCA", treatment_naive = TRUE),
                        cores = 6) |> 
  gene_mutation_rates(covariates = ces.refset.hg38$covariates$ESCA) |> 
  ces_variant(cores = 6)

cesa$selection$selection.1 |>  
  arrange(desc(selection_intensity))

save_cesa(cesa = cesa, file = "output_data/cesa_files/ESCA_cesa.rds")




# CESC ---- 

cesc_maf <- preload_maf(maf = "output_data/CESC.maf",refset = "ces.refset.hg38")

cesa <- CESAnalysis(refset = "ces.refset.hg38") |> 
  load_maf(maf = cesc_maf) |> 
  trinuc_mutation_rates(signature_set = "COSMIC_v3.2",
                        signature_exclusions = 
                          suggest_cosmic_signature_exclusions(cancer_type = "CESC", treatment_naive = TRUE),
                        cores = 6) |> 
  gene_mutation_rates(covariates = ces.refset.hg38$covariates$CESC) |> 
  ces_variant(cores = 6)

cesa$selection$selection.1 |>  
  arrange(desc(selection_intensity))

save_cesa(cesa = cesa, file = "output_data/cesa_files/CESC_cesa.rds")



# BLCA ---- 



blca_maf <- preload_maf(maf = "output_data/BLCA.maf",refset = "ces.refset.hg38")

cesa <- CESAnalysis(refset = "ces.refset.hg38") |> 
  load_maf(maf = blca_maf) |> 
  trinuc_mutation_rates(signature_set = "COSMIC_v3.2",
                        signature_exclusions = 
                          suggest_cosmic_signature_exclusions(cancer_type = "BLCA", treatment_naive = TRUE),
                        cores = 6) |> 
  gene_mutation_rates(covariates = ces.refset.hg38$covariates$bladder) |> 
  ces_variant(cores = 6)

cesa$selection$selection.1 |>  
  arrange(desc(selection_intensity))

save_cesa(cesa = cesa, file = "output_data/cesa_files/BLCA_cesa.rds")


# lihc ---- 

lihc_maf <- preload_maf(maf = "output_data/LIHC.maf",refset = "ces.refset.hg38")

cesa <- CESAnalysis(refset = "ces.refset.hg38") |> 
  load_maf(maf = lihc_maf) |> 
  trinuc_mutation_rates(signature_set = "COSMIC_v3.2",
                        signature_exclusions = 
                          suggest_cosmic_signature_exclusions(cancer_type = "LIHC", treatment_naive = TRUE),
                        cores = 6) |> 
  gene_mutation_rates(covariates = ces.refset.hg38$covariates$LIHC) 


# changed to be all, not just recurrent 
cesa <- cesa |> ces_variant(cores = 6)

cesa$selection$selection.1 |>  
  arrange(desc(selection_intensity))


save_cesa(cesa = cesa, file = "output_data/cesa_files/LIHC_cesa.rds")

