# Description: get a gene_table for coloc using the bonferroni or FDR correction 
            # and only keep genes with significant correlation for a neurodegen trait

# Packages -------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(data.table)

# Arguments ------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

bonf_sign_results <- read.table(here(project_dir, "lava_results", "tsv_univ_bivar_QTLs", "bonferroni_significant_onek1k_bivar_window_100000_20220711.tsv"), 
                           sep = "\t", header = T) %>%
  mutate(trait = stringr::str_remove(phen1, "_[:alnum:]+"),
         cell_type = stringr::str_remove(phen2, "_onek1k_[:alnum:]+")) %>%
  dplyr::select(chr,
                gene = gene_name,
                trait,
                cell_type,
                start,
                end = stop,
                ensembl_id = gene_locus)

fdr_sign_results <- read.table(here(project_dir, "lava_results", "tsv_univ_bivar_QTLs", "FDR_significant_onek1k_bivar_window_100000_20220711.tsv"), 
                                sep = "\t", header = T) %>%
  mutate(trait = stringr::str_remove(phen1, "_[:alnum:]+"),
         cell_type = stringr::str_remove(phen2, "_onek1k_[:alnum:]+")) %>%
  dplyr::select(chr,
                gene = gene_name,
                trait,
                cell_type,
                start,
                end = stop,
                ensembl_id = gene_locus)

neurodegen_genes <- bonf_sign_results %>%
  filter(., trait %in% c("ALS", "AD", "PD", "LBD", "FTD", "MS")) %>%
  select(ensembl_id)

neurodegen_genes_fdr <- fdr_sign_results %>%
  filter(., trait %in% c("ALS", "AD", "PD", "LBD", "FTD", "MS")) %>%
  select(ensembl_id)

bonf_sign_results <- bonf_sign_results %>%
  filter(., ensembl_id %in% neurodegen_genes$ensembl_id)

fdr_sign_results <- fdr_sign_results %>%
  filter(., ensembl_id %in% neurodegen_genes_fdr$ensembl_id)

write.table(bonf_sign_results, here(project_dir, "colocalization", "gene_table_coloc_bonferroni.txt"), sep = "\t", row.names = F, quote = F)
write.table(fdr_sign_results, here(project_dir, "colocalization", "gene_table_coloc_FDR.txt"), sep = "\t", row.names = F, quote = F)

