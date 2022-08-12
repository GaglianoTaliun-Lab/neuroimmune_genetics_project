# Get file info file for LAVA  - for GWAS traits and QTLs combined

#------------------Packages
library(dplyr)
library(here)
library(stringr)

#------------------Arguments
project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# read sample sizes for each GWAS
N_samples <- read.table(here(project_dir,"ALL_case_control_N.csv"), sep=",", header = T) %>% 
  arrange(., trait, by_group=F) %>%
  dplyr::select(
    phenotype = trait,
    cases = n_cases,
    controls = n_controls)

# get path for GWAS lava sumstats
file_path_GWAS <- data.frame(filename = list.files(
    path = here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS"),
    full.names = T,
    recursive = T,
    pattern = ".lava.gz"
  )
)

# filter by retaning only those phenotypes in the latest LAVA run for GWAS and  merge both
input_info_GWAS <- file_path_GWAS %>% 
  mutate(phenotype = 
           stringr::str_remove(filename, "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project/GWAS_summary_statistics/LAVA_sumstats/GWAS/"),
         phenotype =
           stringr::str_remove(phenotype, ".lava.gz")
  ) %>%
  inner_join(., N_samples, by = "phenotype") %>%
  select(.,
         phenotype,
         cases,
         controls,
         filename)

# get path for QTL lava sumstats
file_path_QTL <- data.frame(filename = list.files(
    path = here(project_dir, "GWAS_summary_statistics", "LAVA_sumstats", "1MscBloodNL"),
    full.names = T,
    recursive = T,
    pattern = ".lava.gz"
  )
)

# extract QTL names from list of files
QTL_names <- file_path_QTL %>%
	       mutate(
		 filename = 
                   stringr::str_remove(filename, here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","1MscBloodNL","")),
	         filename = 
		   stringr::str_remove(filename, ".lava.gz")
               )

# fill in data frame with NAs in n_cases and n_controls for QTLs
QTL_names <- data.frame(phenotype = QTL_names$filename, cases = NA, controls = NA)

# merge names and file paths
input_info_QTL <- cbind(QTL_names, file_path_QTL)

# remove "_gene_filtered" from QTL names for output name
# out_qtl <- QTL_names$phenotype

#-----------------Merge and save table

# for (i in 1:nrow(input_info_QTL)){

#  rbind(input_info_GWAS, input_info_QTL[i,], row.names = NULL) %>%
#  write.table(
#    .,
#    here(project_dir,"info_files",stringr::str_c("input.info.gene_gwas.",out_qtl[i],".txt")),
#    sep = "\t", row.names = F, quote = F
#  )

# }

rbind(input_info_GWAS, input_info_QTL, row.names = NULL) %>%
  write.table(
    .,
    here(project_dir,"info_files",stringr::str_c("input.info.gene_gwas.txt")),
    sep = "\t", row.names = F, quote = F
  )
