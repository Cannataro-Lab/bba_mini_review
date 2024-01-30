# download data

library(cancereffectsizeR)


cancereffectsizeR::get_TCGA_project_MAF(project = "LUSC", filename = "output_data/LUSC.maf")
cancereffectsizeR::get_TCGA_project_MAF(project = "ESCA", filename = "output_data/ESCA.maf")
cancereffectsizeR::get_TCGA_project_MAF(project = "CESC", filename = "output_data/CESC.maf")
cancereffectsizeR::get_TCGA_project_MAF(project = "BLCA", filename = "output_data/BLCA.maf")
cancereffectsizeR::get_TCGA_project_MAF(project = "THCA",filename = "output_data/THCA.maf")
cancereffectsizeR::get_TCGA_project_MAF(project = "LUAD", filename = "output_data/LUAD.maf")


