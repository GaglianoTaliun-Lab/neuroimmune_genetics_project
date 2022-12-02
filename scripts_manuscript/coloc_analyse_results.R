# Packages -------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

# Arguments ----------------------------------------------------------------------

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project/colocalization"

significant_results <- list()

# which results? FDR or bonferroni?
res_type = "FDR"

# Read files ---------------------------------------------------------------------

coloc_summ_files <- list.files(here(project_dir, stringr::str_c("coloc_results_",res_type), "summary_all"), full.names = F) %>%
  stringr::str_remove(., "coloc_summary_")

coloc_summ_path <- list.files(here(project_dir, stringr::str_c("coloc_results_", res_type), "summary_all"), full.names = T)

gene_table <- read.table(here(project_dir, stringr::str_c("gene_table_coloc_", res_type,".txt")), sep = "\t", header = T)

# Main ---------------------------------------------------------------------------

# loop across each summary file
for (i in 1:length(coloc_summ_path)) {
  
  # read the summary results
  result <- read.table(coloc_summ_path[i], sep = "\t", header = F) %>%
    rename(header = V1, value = V2)
  
  # save summary table in new folder if H4 >= 0.8
  if (result[6,2] >= 0.8) {
    
    write.table(result, here(project_dir, stringr::str_c("coloc_results_",res_type), "summary_H4_0.8", coloc_summ_files[[i]]), 
                sep = "\t", quote = F, row.names = F)
    
    # read results table for those tests with H4 >= 0.8 and filter to keep only SNPs with a high posterior probability (>= 0.7)
    sign_snps <- read.table(here(project_dir, stringr::str_c("coloc_results_",res_type), "results_all", str_c("coloc_results_", coloc_summ_files[[i]])),
                          sep = "\t", header = T) %>%
      filter(., SNP.PP.H4 >= 0.7)
    
    # save filtered results table in new folder
    if (nrow(sign_snps > 0)) {
      
      write.table(sign_snps, here(project_dir, stringr::str_c("coloc_results_",res_type), "results_H4_0.8_PP_SNP_0.7", coloc_summ_files[[i]]),
                  sep = "\t", quote = F, row.names = F)
      
    }
    
  }
  
  # save summary table in new folder if H3 >= 0.8
  if (result[5,2] >= 0.8) {
    
    write.table(result, here(project_dir, stringr::str_c("coloc_results_",res_type), "summary_H3_0.8", coloc_summ_files[[i]]), 
                sep = "\t", quote = F, row.names = F)
    
  }
  
}

# create a table of significant results (H4 >= 0.8) -----------------------------------------------------------------------------

coloc_sign_files <- as.data.frame(list.files(here(project_dir, stringr::str_c("coloc_results_",res_type), "summary_H4_0.8"), full.names = F))
colnames(coloc_sign_files) <- "files"
coloc_sign_files <- coloc_sign_files %>%
  mutate(files = stringr::str_remove(files, ".tsv")) %>%
  separate(files, c("trait", "cell_type", NA, "ensembl_id"), sep = "_")

sign_coloc <- left_join(coloc_sign_files, gene_table, by = c("trait", "cell_type", "ensembl_id")) %>%
  mutate(., nsnps = NA, H4 = NA)

for (i in 1:nrow(sign_coloc)) {
    
  res_sign_coloc <- read.table(here(project_dir, stringr::str_c("coloc_results_",res_type), "summary_H4_0.8", stringr::str_c(sign_coloc$trait[i], "_", sign_coloc$cell_type[i], "_onek1k_", sign_coloc$ensembl_id[i], ".tsv")),
                               sep = "\t", header = T)
  
  sign_coloc$H4[i] <- res_sign_coloc[6,2]
  sign_coloc$nsnps[i] <- res_sign_coloc[1,2]
  
}

write.table(sign_coloc, here(project_dir, stringr::str_c("coloc_results_",res_type), "summary_H4_0.8", "summary_of_results_H4_0.8_v2.tsv"), sep = "\t", row.names = F, quote = F)

# create a table of H3 significant results (H3 >= 0.8) -----------------------------------------------------------------------------

coloc_sign_files <- as.data.frame(list.files(here(project_dir, stringr::str_c("coloc_results_",res_type), "summary_H3_0.8"), full.names = F))
colnames(coloc_sign_files) <- "files"
coloc_sign_files <- coloc_sign_files %>%
  mutate(files = stringr::str_remove(files, ".tsv")) %>%
  separate(files, c("trait", "cell_type", NA, "ensembl_id"), sep = "_")

sign_coloc <- left_join(coloc_sign_files, gene_table, by = c("trait", "cell_type", "ensembl_id")) %>%
  mutate(., nsnps = NA, H3 = NA)

for (i in 1:nrow(sign_coloc)) {
  
  res_sign_coloc <- read.table(here(project_dir, stringr::str_c("coloc_results_",res_type), "summary_H3_0.8", stringr::str_c(sign_coloc$trait[i], "_", sign_coloc$cell_type[i], "_onek1k_", sign_coloc$ensembl_id[i], ".tsv")),
                               sep = "\t", header = T)
  
  sign_coloc$H3[i] <- res_sign_coloc[5,2]
  sign_coloc$nsnps[i] <- res_sign_coloc[1,2]
  
}

write.table(sign_coloc, here(project_dir, stringr::str_c("coloc_results_",res_type), "summary_H3_0.8", "summary_of_results_H3_0.8.tsv"), sep = "\t", row.names = F, quote = F)
