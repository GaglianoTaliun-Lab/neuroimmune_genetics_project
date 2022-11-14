# Description: obtain list of test loci to replicate the AD LAVA with onek1k results and pQTLs on Wightman et al. and Kunkle et al.

# Load packages -----------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(data.table)

# Set arguments -----------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# Load data ---------------------------------------------------------------

FDR_significant_onek1k <- read.table(here(project_dir, "lava_results", "tsv_univ_bivar_QTLs", "FDR_significant_onek1k_bivar_window_100000_20220711.tsv"),
  sep = "\t", header = T)

FDR_significant_pqtl <- read.table(here(project_dir, "lava_results", "tsv_univ_bivar_QTLs", "FDR_significant_pQTLs_bivar_window_100000_20220803.tsv"),
  sep = "\t", header = T)

onek1k_test_loci <- read.table(here(project_dir, "test_loci", "window_100000", "onek1k_filtered.loci"), 
  sep = "\t", header = T)

pqtl_test_loci <- read.table(here(project_dir, "test_loci", "window_100000", "pqtls_filtered.loci"),
  sep = "\t", header = T)

# Wrangle data ---------------------------------------------------------------

genes_AD_onek1k <- FDR_significant_onek1k %>%
  filter(., phen1 == "AD_schwartzentruber2021") %>%
  dplyr::select(gene_locus)

AD_replicate_onek1k <- onek1k_test_loci %>%
  filter(., LOC %in% genes_AD_onek1k$gene_locus)

genes_AD_pqtl <- FDR_significant_pqtl %>%
  filter(., phen1 == "AD_schwartzentruber2021") %>%
  dplyr::select(gene_locus)

AD_replicate_pqtl <- pqtl_test_loci %>%
  filter(., LOC %in% genes_AD_pqtl$gene_locus)

# Write new test loci file ---------------------------------------------------------------

write.table(AD_replicate_onek1k, here(project_dir, "test_loci", "window_100000", "AD_replication_onek1k.loci"),
  sep = "\t", row.names = F, col.names = T, quote = F) 

write.table(AD_replicate_pqtl, here(project_dir, "test_loci", "window_100000", "AD_replication_pqtls.loci"),
  sep = "\t", row.names = F, col.names = T, quote = F)
