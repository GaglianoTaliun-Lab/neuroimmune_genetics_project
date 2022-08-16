# Description: run colocalization analysis between QTLs and GWAS sumstats after LAVA

# Packages -------------------------------------------------------

library(coloc)
library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(colochelpR)

# Arguments ----------------------------------------------------------------------

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# default options:
p1 = 1e-04
p2 = 1e-04
p12 = 1e-05

# Read files -------------------------------------------------------

gene_table <- read.table(here(project_dir, "colocalization", "gene_table_coloc.txt"),
                         sep = "\t", header = T)

GWAS_for_coloc <- list.files(here(project_dir, "colocalization", "GWAS_regional_sumstats"), pattern = "*.tsv", all.files = T, full.names = T)

QTLs_for_coloc <- list.files(here(project_dir, "colocalization", "QTL_regional_sumstats"), pattern = "*.tsv", all.files = T, full.names = T)

# Run coloc ----------------------------------------------------------------------------------

coloc_results_summ <- list()
coloc_results_res <- list()
results_names <- array()

for (i in 1:nrow(gene_table)){
  
  GWAS_file_path = stringr::str_c(here(project_dir, "colocalization", "GWAS_regional_sumstats"), "/", gene_table$trait[i], "_", gene_table$gene[i], ".tsv")
  QTL_file_path = stringr::str_c(here(project_dir, "colocalization", "QTL_regional_sumstats"), "/onek1k_", gene_table$cell_type[i], "_", gene_table$gene[i], ".tsv")
  
  GWAS_df <- read.table(GWAS_file_path, sep = "\t", header= T) %>%
    filter(., MAF > 0) %>%
    distinct(., SNP, .keep_all = TRUE)
  
  QTL_df <- read.table(QTL_file_path, sep = "\t", header = T) %>%
    filter(., MAF > 0) %>%
    distinct(., SNP, .keep_all = TRUE)

  if ("beta" %in% names(GWAS_df)) {
    
    coloc_results <- coloc.abf(dataset1 = list(type = "quant",
                                               snp = GWAS_df$SNP,
                                               beta = GWAS_df$beta,
                                               varbeta = GWAS_df$varbeta,
                                               pvalues = GWAS_df$pvalues,
                                               MAF = GWAS_df$MAF,
                                               N = GWAS_df$N),
                               dataset2 = list(type = "quant",
                                               snp = QTL_df$SNP,
                                               beta = QTL_df$beta,
                                               varbeta = QTL_df$varbeta,
                                               pvalues = QTL_df$pvalues,
                                               MAF = QTL_df$MAF,
                                               N = QTL_df$N),
                               p1 = p1, p2 = p2, p12 = p12
    )
    
    coloc_results_summ[[i]] <- coloc_results$summary
    coloc_results_res[[i]] <- coloc_results$results
    results_names[i] <- str_c(GWAS_df$GWAS[1],"_",QTL_df$eQTL_dataset[1])
    
  } else if ("OR" %in% names(GWAS_df)) {
    
    coloc_results <- coloc.abf(dataset1 = list(type = "quant",
                                               snp = GWAS_df$SNP,
                                               OR = GWAS_df$OR,
                                               pvalues = GWAS_df$pvalues,
                                               MAF = GWAS_df$MAF,
                                               N = GWAS_df$N),
                               dataset2 = list(type = "quant",
                                               snp = QTL_df$SNP,
                                               beta = QTL_df$beta,
                                               varbeta = QTL_df$varbeta,
                                               pvalues = QTL_df$pvalues,
                                               MAF = QTL_df$MAF,
                                               N = QTL_df$N),
                               p1 = p1, p2 = p2, p12 = p12
    )
    
    coloc_results_summ[[i]] <- coloc_results$summary
    coloc_results_res[[i]] <- coloc_results$results
    results_names[i] <- str_c(GWAS_df$GWAS[1],"_",QTL_df$eQTL_dataset[1])
    
  } else next
  
}

for (i in 1:length(results_names)) {
  
  summ_out <- as.data.frame(coloc_results_summ[[i]]) %>% drop_na()
  res_out <- as.data.frame(coloc_results_res[[i]]) %>% drop_na()
  write.table(summ_out, here(project_dir, "colocalization", "coloc_results", str_c("coloc_summary_", results_names[i],".tsv")), sep = "\t", row.names = T, quote = F, col.names = F)
  write.table(res_out, here(project_dir, "colocalization", "coloc_results", str_c("coloc_results_", results_names[i],".tsv")), sep = "\t", row.names = F, quote = F)
                              
}


