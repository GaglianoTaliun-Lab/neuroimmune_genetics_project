# Description: script to retrieve information about genome-wide significant 
              # signals in GWAS from significant regional correlations

# Packages -------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(data.table)

# Arguments ------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# Load files -----------------------------------------------------

sign_pqtls_fdr <- read.table(here(project_dir, "lava_results", "tsv_univ_bivar_QTLs", "FDR_significant_pQTLs_bivar_window_100000_20220803.tsv"), sep = "\t", header = T)
sign_onek1k_fdr <- read.table(here(project_dir, "lava_results", "tsv_univ_bivar_QTLs", "FDR_significant_onek1k_bivar_window_100000_20220711.tsv"), sep = "\t", header = T)

path_GWAS_sumstats <- list.files(path = here(project_dir, "GWAS_summary_statistics", "LAVA_sumstats", "GWAS"), pattern = "*.lava.gz", full.names = T)

# Main -----------------------------------------------------------

# create new column
sign_pqtls_fdr <- sign_pqtls_fdr %>%
  mutate(GWAS_significant = NA)

sign_onek1k_fdr <- sign_onek1k_fdr %>%
  mutate(GWAS_significant = NA)

# evaluate each line from the sign_pqtls_fdr 
for (i in 1:nrow(sign_pqtls_fdr)) {

  cat("---- Evaluating line", i, "out of", nrow(sign_pqtls_fdr), ".\n")

  GWAS <- sign_pqtls_fdr$phen1[i]
  start <- sign_pqtls_fdr$start[i]
  stop <- sign_pqtls_fdr$stop[i]
  
  GWAS_sumstats <- fread(here(project_dir, "GWAS_summary_statistics", "LAVA_sumstats", "GWAS", stringr::str_c(GWAS,".lava.gz"))) %>%
    filter(., BP >= start & BP <= stop) %>%
    filter(., P < 5e-08)
  
  if (nrow(GWAS_sumstats) > 0) {
    sign_pqtls_fdr$GWAS_significant[i] <- "YES"
  } else if (nrow(GWAS_sumstats) == 0) {
    sign_pqtls_fdr$GWAS_significant[i] <- "NO"
  }
  
}

cat("---- Done evaluating the significant signals for pQTLs!\n")

write.table(sign_pqtls_fdr, here(project_dir, "lava_results", "tsv_univ_bivar_QTLs", "FDR_significant_pQTLs_bivar_window_100000_GWASsignificant_20221001.tsv"), sep = "\t", quote = F, row.names = F)

# evaluate each line from the sign_onek1k_fdr 
for (i in 1:nrow(sign_onek1k_fdr)) {

  cat("---- Evaluating line", i, "out of", nrow(sign_onek1k_fdr), ".\n")

  GWAS <- sign_onek1k_fdr$phen1[i]
  start <- sign_onek1k_fdr$start[i]
  stop <- sign_onek1k_fdr$stop[i]
  
  GWAS_sumstats <- fread(here(project_dir, "GWAS_summary_statistics", "LAVA_sumstats", "GWAS", stringr::str_c(GWAS,".lava.gz"))) %>%
    filter(., BP >= start & BP <= stop) %>%
    filter(., P < 5e-08)

  if (nrow(GWAS_sumstats) > 0) {
    sign_onek1k_fdr$GWAS_significant[i] <- "YES"
  } else if (nrow(GWAS_sumstats) == 0) {
    sign_onek1k_fdr$GWAS_significant[i] <- "NO"
  }  
  
}

cat("---- Done evaluating the significant signals for eQTLs!\n") 

write.table(sign_onek1k_fdr, here(project_dir, "lava_results", "tsv_univ_bivar_QTLs", "FDR_significant_onek1k_bivar_window_100000_GWASsignificant_20221001.tsv"), sep = "\t", quote = F, row.names = F)


