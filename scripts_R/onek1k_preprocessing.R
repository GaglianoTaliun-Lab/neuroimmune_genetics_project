# Description: load and wrangle Onek1k beta summary statistics for LAVA per gene

# Load packages -----------------------------------------------------------

library(here)
library(data.table)
library(stringr)
library(tidyverse)
library(rutils)
library(R.utils)
library(utils)

# Arguments ---------------------------------------------------------------

args = commandArgs(TRUE)

cell_type_file_name <- args[1]

cell_type_in_file <- args[2]

cat("Argument 1 = ", cell_type_file_name,". Argument 2 = ", cell_type_in_file, ".\n")

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

write.file.function <- function(df_out) {
  
  if (!is.na(nrow(df_out))) {
    
    out_name=here(project_dir, "GWAS_summary_statistics", "LAVA_sumstats", "onek1k", paste0(cell_type_file_name,"_onek1k_",df_out$GENE[1],".lava.gz"))
    
    if (!is.na(out_name)) {
      fwrite(df_out, file = out_name, sep = "\t")
      
    }
    
  }

}

# List all files ---------------------------------------------------------------

# read rho file sumstats to add alleles
rho_file <- read.table(here(project_dir,"GWAS_summary_statistics","onek1k","onek1k_eqtl_dataset.tsv"), sep = "\t", header = TRUE) %>%
  filter(., CELL_ID == cell_type_in_file) %>%
  select(SNP = RSID,
         CHR,
         POS,
         A1,
         A2,
         gene = GENE,
         GENE_ID)

# read all chromosomes from specific cell type  
all_files_path <- list.files(path = here(project_dir, "GWAS_summary_statistics", "onek1k", "OneK1K_matrix_eQTL_results"), 
           pattern = str_c("^",cell_type_file_name), full.names = TRUE, all.files = TRUE)

all_files <- list.files(path = here(project_dir, "GWAS_summary_statistics", "onek1k", "OneK1K_matrix_eQTL_results"), 
                             pattern = str_c("^",cell_type_file_name), full.names = FALSE, all.files = TRUE)

# Main loop
for (i in 1:length(all_files)){
  
  one_file <- read.table(all_files_path[i], sep = "\t", header = TRUE) %>%
    rename(., pvalue = p.value) %>%
    rename(., tstat = t.stat)
  
  # get list of unique genes
  unique_genes <- one_file %>%
    filter(., pvalue < 5e-08) %>%
    dplyr::select(gene) %>%
    distinct(., gene)

  cat("For ", all_files[i], " there are ", nrow(unique_genes), " unique genes.\n")
  
  # prepare lava.gz file with only significant genes  
  one_file %>%
    dplyr::filter(gene %in% unique_genes$gene) %>%
    mutate(SE = beta/tstat, N = 982) %>%
    left_join(., rho_file, by = c("SNP", "gene")) %>%
    dplyr::select(.,
                  SNP,
                  GENE = GENE_ID,
                  CHR,
                  BP = POS,
                  A1 = A2,
                  A2 = A1, # A2 is the effect allele (confirmed by email from authors)
                  BETA = beta,
                  SE,
                  P = pvalue,
                  N) %>%
    filter(., !is.na(GENE)) %>%
    group_split(GENE) %>%
    as.list(.) %>%

    lapply(., write.file.function)

}
  
