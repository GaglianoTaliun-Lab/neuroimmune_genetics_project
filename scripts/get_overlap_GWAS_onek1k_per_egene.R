# Description: get sample overlaps for eQTL/GWAS traits

# Load packages -----------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(qdapTools)

# Set arguments -----------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

out_dir <- here(project_dir, "sample_overlap")

args <- commandArgs(TRUE)
cell_type = args[1]

gwas <- c("AD_schwartzentruber2021", "LBD_chia2021", "PD_nalls2019", "ALS_vanrheenen2021",
          "FTD_ferrari2014", "MS_imsgc2019", "CD_delange2017", "UC_delange2017", "SCZ_pardinas2018") %>% sort()

eqtl_gene <-
  list.files(
    path = here(project_dir, "GWAS_summary_statistics", "LAVA_sumstats","onek1k"),
    pattern = "lava.gz"
  ) %>%
  stringr::str_subset(., str_c("^",cell_type)) %>%
  basename() %>%
  str_remove(".lava.gz")

cat("The number of eQTL genes for", cell_type, "is:", length(eqtl_gene), ".\n")

phenotypes <- c(gwas, eqtl_gene)

# Load ldsc results -------------------------------------------------------

gwas_rg <-
  read.table(
    file = here(project_dir, "ldsc_corr", "ldsc_correlations.txt"),
    sep = "\t",
    header = T
  )

# Main --------------------------------------------------------------------

# Generate all possible permutations with repetitions
permutations <-
  tibble(
    p1 = c(gwas, eqtl_gene),
    p2 = c(gwas, eqtl_gene)
  ) %>%
  tidyr::expand(p1, p2)

# LDSC already run for gwas phenotypes
# Just need to add 0 in where eqtl phenotypes are used
# We are assuming that no sample overlaps exist between gwas and eqtl cohorts
all_rg <-
  permutations %>%
  dplyr::left_join(
    gwas_rg
  ) %>%
  dplyr::select(p1, p2, gcov_int) %>%
  dplyr::mutate(
    gcov_int =
      case_when(
        p1 %in% eqtl_gene & !p2 %in% eqtl_gene ~ 0,
        !p1 %in% eqtl_gene & p2 %in% eqtl_gene ~ 0,
        p1 %in% eqtl_gene & p2 %in% eqtl_gene & p1 != p2 ~ 0,
        p1 %in% eqtl_gene & p2 %in% eqtl_gene & p1 == p2 ~ 1,
        TRUE ~ gcov_int
      )
  )

###### Creating sample overlap matrix ######

n <- length(phenotypes)
covar_matrix <- matrix(NA,n,n)
rownames(covar_matrix) <- colnames(covar_matrix) <- phenotypes

for(i in phenotypes) {
  for(j in phenotypes) {
    
    covar_matrix[i,j] <-
      all_rg %>%
      dplyr::filter(p1 == i, p2 == j) %>%
      .[["gcov_int"]]
    
  }
}

# Sometimes there might be small differences in gcov_int depending on which phenotype was analysed as the outcome / predictor
if (!all(t(covar_matrix)==covar_matrix)) {
  covar_matrix[lower.tri(covar_matrix)] <- t(covar_matrix)[lower.tri(covar_matrix)]
}

# Standardise the matrix
covar_matrix <-
  covar_matrix %>%
  cov2cor() %>%
  round(digits = 5)

# Save data ---------------------------------------------------------------

write.table(
  covar_matrix,
  file = file.path(out_dir,
                   str_c(cell_type,"_onek1k_gwas_sample_overlap.txt")),
  quote = F
)
