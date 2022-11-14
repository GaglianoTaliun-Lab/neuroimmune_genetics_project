# Description: subset of the script 'get_sample_overlaps.R' from Regina,
# to process after running 'run_ldsc.sh' across all pairs of traits and obtain a correlation matrix
# and the sample overlaps.

# Load packages -----------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(data.table)
library(readr)

# Set arguments -----------------------------------------------------------

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# phenotypes to include in the matrix
phenotypes_prox <- read.table(here(project_dir, "GWAS_summary_statistics", "list_of_GWAS_studies.txt"), sep = "\t", header = F) %>% .[,1] %>% as.array() %>% sort()
phenotypes_noprox <- read.table(here(project_dir, "GWAS_summary_statistics", "list_of_GWAS_studies_no_proxies.txt"), sep = "\t", header = F) %>% .[,1] %>% as.array() %>% sort()

# read and write input/output directories:
gwas_dir <- here(project_dir, "GWAS_summary_statistics", "ldsc_sumstats")

out_dir <- here(project_dir,"ldsc_corr")
out_dir2 <- here(project_dir,"sample_overlap")

# number of lines to skip in the ldsc output:
n_skip=60
# number of lines to read from the output:
n_read=1

###### Extracting LDSC results ######

# Extract all outputs into single file
file_paths <-
  list.files(
    out_dir,
    pattern = "_rg.log",
    full.names = T
  )

files <-
  vector(mode = "list",
         length = length(file_paths))

for(i in 1:length(file_paths)){
  
  files[[i]] <-
    read_table(
      file = file_paths[i],
      skip = n_skip, # Number of lines to skip in log file
      n_max = n_read # Number of lines to read
    )
  
}

# extract rg values in matrix
all_rg <-
  files %>%
  rbindlist(., fill = TRUE) %>%
  dplyr::mutate(
    p1 = basename(p1) %>%
      stringr::str_remove(".sumstats.gz") %>%
      stringr::str_remove(str_c(gwas_dir,"/")),
    p2 = basename(p2) %>%
      stringr::str_remove(".sumstats.gz") %>%
      stringr::str_remove(str_c(gwas_dir,"/"))
  ) %>%
  select(
    p1,
    p2,
    rg,
    se,
    z,
    p,
    h2_obs,
    h2_obs_se,
    h2_int,
    h2_int_se,
    gcov_int,
    gcov_int_se
  )

###### Creating sample overlap matrix by extracting the intercept from LDSC results ######

n <- length(phenotypes_prox)
n_noprox <- length(phenotypes_noprox)
covar_matrix <- list(matrix(NA,n,n), matrix(NA,n_noprox,n_noprox))
phenotypes <- list(phenotypes_prox, phenotypes_noprox)

for (k in 1:length(phenotypes)){
  
  rownames(covar_matrix[[k]]) <- colnames(covar_matrix[[k]]) <- phenotypes[[k]]
  
  for(i in phenotypes[[k]]) {
    for(j in phenotypes[[k]]) {
      
      covar_matrix[[k]][i,j] <-
        all_rg %>%
        dplyr::filter(p1 == i, p2 == j) %>%
        .[["gcov_int"]]
      
    }
  }
  
  # Sometimes there might be small differences in gcov_int depending on which phenotype was analysed as the outcome / predictor
  if (!all(t(covar_matrix[[k]])==covar_matrix[[k]])) {
    covar_matrix[[k]][lower.tri(covar_matrix[[k]])] <- t(covar_matrix[[k]])[lower.tri(covar_matrix[[k]])]
  }
  
  # Standardise the matrix
  covar_matrix[[k]] <-
    covar_matrix[[k]] %>%
    cov2cor() %>%
    round(digits = 5)
  
}

# Save data ---------------------------------------------------------------

write.table(
  all_rg,
  file = here(out_dir, stringr::str_c("ldsc_correlations.txt")),
  sep = "\t",
  row.names = F,
  quote = F
)

write.table(
  covar_matrix[[1]],
  file = here(out_dir2, stringr::str_c("sample_overlap.txt")),
  quote = F,
  row.names = phenotypes_prox,
  sep = "\t"
)

write.table(
  covar_matrix[[2]],
  file = here(out_dir2, stringr::str_c("sample_overlap_no_proxies.txt")),
  quote = F,
  row.names = phenotypes_noprox,
  sep = "\t"
)

