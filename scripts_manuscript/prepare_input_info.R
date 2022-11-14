# Get file info file for LAVA  - only for GWAS traits

library(dplyr)
library(here)

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

N_samples <- read.table(here(project_dir,"ALL_case_control_N.csv"), sep=",", header = T) %>% 
  arrange(., trait, by_group=F) %>%
  dplyr::select(
    phenotype = trait,
    cases = n_cases,
    controls = n_controls)

N_samples_noprox <- read.table(here(project_dir,"ALL_case_control_N_no_proxies.csv"), sep=",", header = T) %>%
  arrange(., trait, by_group=F) %>%
  dplyr::select(
    phenotype = trait,
    cases = n_cases,
    controls = n_controls)


file_path <- data.frame(filename = list.files(
    path = here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS"),
    full.names = T,
    recursive = T,
    pattern = ".lava.gz"
  )
)

input_info <- cbind(N_samples, file_path)
input_info_noprox <- cbind(N_samples_noprox, file_path)

write.table(
  input_info,
  here(project_dir,"info_files","input.info.txt"),
  sep = "\t", row.names = F, quote = F
)

write.table(
  input_info_noprox,
  here(project_dir,"info_files","input.info.no_proxies.txt"),
  sep = "\t", row.names = F, quote = F
)
