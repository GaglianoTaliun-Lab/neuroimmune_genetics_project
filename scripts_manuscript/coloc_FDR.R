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

gene_table <- read.table(here(project_dir, "colocalization", "gene_table_coloc_FDR.txt"),
                         sep = "\t", header = T)

GWAS_for_coloc <- list.files(here(project_dir, "colocalization", "GWAS_regional_sumstats"), pattern = "*.tsv", all.files = T, full.names = T)

QTLs_for_coloc <- list.files(here(project_dir, "colocalization", "QTL_regional_sumstats"), pattern = "*.tsv", all.files = T, full.names = T)

# table to get case/control fraction:
case_control_N <- read.table(here(project_dir, "ALL_case_control_N.csv"), sep = ",", header = T) %>%
  mutate(fraction_cases = n_cases/(n_cases + n_controls),
         disease = str_remove(trait, "_[:alnum:]+")) 

head(case_control_N)

# Run coloc ----------------------------------------------------------------------------------

coloc_results_summ <- list()
coloc_results_res <- list()
results_names <- array()

for (i in 1:nrow(gene_table)){
  
  GWAS_file_path = stringr::str_c(here(project_dir, "colocalization", "GWAS_regional_sumstats"), "/", gene_table$trait[i], "_", gene_table$gene[i], ".tsv")
  QTL_file_path = stringr::str_c(here(project_dir, "colocalization", "QTL_regional_sumstats"), "/coloc_onek1k_", gene_table$cell_type[i], "_", gene_table$gene[i], ".tsv")
  
  # get s (fraction of cases):
  fraction_cases <- filter(case_control_N, disease == gene_table$trait[i]) %>%
    .[1,4]
  
  GWAS_df <- read.table(GWAS_file_path, sep = "\t", header= T) %>%
    filter(., MAF > 0 & MAF < 1) %>%
    distinct(., SNP, .keep_all = TRUE)
  
  QTL_df <- try(read.table(QTL_file_path, sep = "\t", header = T), silent = TRUE)
  
  if (class(QTL_df) != "try-error") {

    if (nrow(QTL_df) > 0) {
    
      QTL_df <- QTL_df %>% filter(., MAF > 0 & MAF < 1) %>%
        distinct(., SNP, .keep_all = TRUE)
    
      # get unique list of SNPs across both datasets
      unique_snps <- inner_join(GWAS_df, QTL_df, by = "SNP") %>% .$SNP 
    
      GWAS_df <- GWAS_df %>%
        filter(., SNP %in% unique_snps)
    
      QTL_df <- QTL_df %>%
        filter(., SNP %in% unique_snps)
    
      if ("beta" %in% names(GWAS_df)) {
      
        coloc_results <- coloc.abf(dataset1 = list(type = "cc",
                                                 snp = GWAS_df$SNP,
                                                 beta = GWAS_df$beta,
                                                 varbeta = GWAS_df$varbeta,
                                                 pvalues = GWAS_df$pvalues,
                                                 MAF = GWAS_df$MAF,
                                                 N = GWAS_df$N,
                                                 s = fraction_cases),
                                   dataset2 = list(type = "quant",
                                                 snp = QTL_df$SNP,
                                                 pvalues = QTL_df$pvalues,
                                                 MAF = QTL_df$MAF,
                                                 N = QTL_df$N),
                                   p1 = p1, p2 = p2, p12 = p12
        )
      
        coloc_results_summ[[i]] <- coloc_results$summary
        coloc_results_res[[i]] <- coloc_results$results
        results_names[i] <- str_c(GWAS_df$GWAS[1],"_",QTL_df$eQTL_dataset[1])
      
        pdf(here(project_dir, "colocalization", "sensitivity_analysis", str_c(results_names[i], ".pdf")))
        cat("Sensitivity analysis for locus <", results_names[i],"> :\t")
        sensitivity(coloc_results,rule="H4 >= 0.8")
        dev.off()
      
        cat("Coloc run #", i, "-- Analysis done!\n")
      
      } else if ("OR" %in% names(GWAS_df)) {
      
        coloc_results <- coloc.abf(dataset1 = list(type = "cc",
                                                 snp = GWAS_df$SNP,
                                                 OR = GWAS_df$OR,
                                                 pvalues = GWAS_df$pvalues,
                                                 MAF = GWAS_df$MAF,
                                                 N = GWAS_df$N,
                                                 s = fraction_cases),
                                   dataset2 = list(type = "quant",
                                                 snp = QTL_df$SNP,
                                                 pvalues = QTL_df$pvalues,
                                                 MAF = QTL_df$MAF,
                                                 N = QTL_df$N),
                                   p1 = p1, p2 = p2, p12 = p12
        )
      
        coloc_results_summ[[i]] <- coloc_results$summary
        coloc_results_res[[i]] <- coloc_results$results
        results_names[i] <- str_c(GWAS_df$GWAS[1],"_",QTL_df$eQTL_dataset[1])
      
        pdf(here(project_dir, "colocalization", "sensitivity_analysis", str_c(results_names[i], ".pdf")))
        cat("Sensitivity analysis for locus <", results_names[i],"> :\t")
        sensitivity(coloc_results,rule="H4 >= 0.8") 
        dev.off()
      
        cat("Coloc run #", i, "-- Analysis done!\n")
      
      } else next
    
    } else next

  } else next
  
}

for (i in 1:length(results_names)) {
  
  summ_out <- as.data.frame(coloc_results_summ[[i]]) %>% drop_na()
  res_out <- as.data.frame(coloc_results_res[[i]]) %>% drop_na()
  
  if (nrow(summ_out) > 0) {
    
    write.table(summ_out, here(project_dir, "colocalization", "coloc_results_FDR", "summary_all", str_c("coloc_summary_", results_names[i],".tsv")), sep = "\t", row.names = T, quote = F, col.names = F)
    write.table(res_out, here(project_dir, "colocalization", "coloc_results_FDR", "results_all", str_c("coloc_results_", results_names[i],".tsv")), sep = "\t", row.names = F, quote = F)
    
  }

}

