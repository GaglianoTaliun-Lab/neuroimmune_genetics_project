# Description: run colocalization analysis with susie between QTLs and GWAS sumstats 
# only for coloc results H4 >= 0.8

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

gene_table <- read.table(here(project_dir, "colocalization", "coloc_results_FDR", "summary_H4_0.8", "summary_of_results_H4_0.8.tsv"),
                         sep = "\t", header = T)

GWAS_for_coloc <- list.files(here(project_dir, "colocalization", "GWAS_regional_sumstats"), pattern = "*.tsv", all.files = T, full.names = T)

QTLs_for_coloc <- list.files(here(project_dir, "colocalization", "QTL_regional_sumstats"), pattern = "*.tsv", all.files = T, full.names = T)

# Run coloc ----------------------------------------------------------------------------------

coloc_results_summ <- list()
coloc_results_res <- list()
results_names <- array()

# for (i in 1:nrow(gene_table)){
for (i in 91:99){

  cat("Starting analysis for row #",i,".\n")
  
  GWAS_file_path = stringr::str_c(here(project_dir, "colocalization", "GWAS_regional_sumstats"), "/", gene_table$trait[i], "_", gene_table$gene[i], ".tsv")
  QTL_file_path = stringr::str_c(here(project_dir, "colocalization", "QTL_regional_sumstats"), "/onek1k_", gene_table$cell_type[i], "_", gene_table$gene[i], ".tsv")
  
  bim_file <- read.table(here(project_dir, "colocalization", "LD_matrices", stringr::str_c("LD_matrix_", gene_table$ensembl_id[i], ".bim")), sep = "\t", header = F)
  colnames(bim_file) <- c("CHR","SNP","V3","BP","REF","ALT")
  cat("First lines of bim file:\n")
  head(bim_file)
  
  LD_matrix <- read.table(here(project_dir, "colocalization", "LD_matrices", stringr::str_c("LD_matrix_", gene_table$ensembl_id[i], ".ld")), sep = "\t", header = F)
  rownames(LD_matrix) <- bim_file[,2]
  colnames(LD_matrix) <- bim_file[,2]
  
  GWAS_df <- read.table(GWAS_file_path, sep = "\t", header= T) %>%
    filter(., MAF > 0) %>%
    distinct(., SNP, .keep_all = TRUE)
  
  if ("beta" %in% names(GWAS_df)) {
    
    GWAS_df <- inner_join(GWAS_df, bim_file, by = c("SNP","CHR","BP")) %>%
      mutate(beta_aligned = 
               case_when(
                 A1 == REF ~ beta * -1,
                 A1 == ALT ~ beta
               )) %>%
      filter(., !is.na(beta_aligned)) %>%
      dplyr::select(GWAS, SNP, CHR, BP, beta = beta_aligned, se, pvalues, A1, A2, N, varbeta, MAF)
    
  } else if ("OR" %in% names(GWAS_df)) {
    
    GWAS_df <- inner_join(GWAS_df, bim_file, by = c("SNP","CHR","BP")) %>%
      mutate(OR_aligned = 
               case_when(
                 A1 == REF ~ OR * -1,
                 A1 == ALT ~ OR
               )) %>%
      filter(., !is.na(OR_aligned)) %>%
      dplyr::select(GWAS, SNP, CHR, BP, OR = OR_aligned, se, pvalues, A1, A2, N, MAF)
    
  } else next

  cat("First lines of GWAS_df:\n")
  head(GWAS_df)
  
  QTL_df <- read.table(QTL_file_path, sep = "\t", header = T) %>%
    filter(., MAF > 0) %>%
    distinct(., SNP, .keep_all = TRUE)
  
  QTL_df <- inner_join(QTL_df, bim_file, by = c("SNP","CHR","BP")) %>%
    mutate(beta_aligned = 
             case_when(
               A1 == REF ~ beta * -1,
               A1 == ALT ~ beta
             )) %>%
    filter(., !is.na(beta_aligned)) %>%
    dplyr::select(eQTL_dataset, gene, SNP, CHR, BP, beta = beta_aligned, se, pvalues, A1, A2, N, varbeta, MAF)
  
  cat("First lines of QTL_df:\n")
  head(QTL_df)
  
  LD_GWAS <- as.matrix(LD_matrix[GWAS_df$SNP, GWAS_df$SNP])
  
  cat("First lines of LD_GWAS:\n")
  # head(LD_GWAS)
  
  LD_QTL <- as.matrix(LD_matrix[QTL_df$SNP, QTL_df$SNP])
  
  cat("First lines of LD_QTL:\n")
  # head(LD_QTL)
  
  if ("beta" %in% names(GWAS_df)) {
    
    dataset1 = list(type = "cc",
                    snp = GWAS_df$SNP,
                    beta = GWAS_df$beta,
                    varbeta = GWAS_df$varbeta,
                    pvalues = GWAS_df$pvalues,
                    MAF = GWAS_df$MAF,
                    N = GWAS_df$N,
                    LD = LD_GWAS)
    
    dataset2 = list(type = "quant",
                    snp = QTL_df$SNP,
                    beta = QTL_df$beta,
                    varbeta = QTL_df$varbeta,
                    pvalues = QTL_df$pvalues,
                    MAF = QTL_df$MAF,
                    N = QTL_df$N,
                    LD = LD_QTL)
    
    S1 <- try(runsusie(dataset1))
    
    if (inherits(S1, "try-error")) {
      next
    }
    
    S2 <- try(runsusie(dataset2))
    
    if (inherits(S1, "try-error")) {
      next
    }    
    coloc_results = coloc.susie(S1, S2)
    
    coloc_results_summ[[i]] <- coloc_results$summary
    coloc_results_res[[i]] <- coloc_results$results
    results_names[i] <- str_c(GWAS_df$GWAS[1],"_",QTL_df$eQTL_dataset[1])
    
  } else if ("OR" %in% names(GWAS_df)) {
    
    dataset1 = list(type = "quant",
                    snp = GWAS_df$SNP,
                    OR = GWAS_df$OR,
                    pvalues = GWAS_df$pvalues,
                    MAF = GWAS_df$MAF,
                    N = GWAS_df$N,
                    LD = LD_GWAS)
    
    dataset2 = list(type = "quant",
                    snp = QTL_df$SNP,
                    beta = QTL_df$beta,
                    varbeta = QTL_df$varbeta,
                    pvalues = QTL_df$pvalues,
                    MAF = QTL_df$MAF,
                    N = QTL_df$N,
                    LD = LD_QTL)
    
    S1 <- try(runsusie(dataset1))
    
    if (inherits(S1, "try-error")) {
      next
    }
    
    S2 <- try(runsusie(dataset2))
    
    if (inherits(S1, "try-error")) {
      next
    }

    coloc_results = coloc.susie(S1, S2)
    
    coloc_results_summ[[i]] <- coloc_results$summary
    coloc_results_res[[i]] <- coloc_results$results
    results_names[i] <- str_c(GWAS_df$GWAS[1],"_",QTL_df$eQTL_dataset[1])
    
  } else next
  
  cat("Finished	analysis for row #",i,".\n")

}

for (i in 1:length(results_names)) {
  
  summ_out <- as.data.frame(coloc_results_summ[[i]]) %>% drop_na()
  res_out <- as.data.frame(coloc_results_res[[i]]) %>% drop_na()

  if (nrow(summ_out) > 0) {
    
    cat("Writing summary table for row #",i,".\n")
    write.table(summ_out, here(project_dir, "colocalization", "coloc_results_FDR", "coloc_with_susie", "summary_all", str_c("coloc.susie_summary_", results_names[i],".tsv")), sep = "\t", row.names = T, quote = F, col.names = F)

  }
  
  if (nrow(res_out) > 0) {

    cat("Writing summary table for row #",i,".\n")
    write.table(res_out, here(project_dir, "colocalization", "coloc_results_FDR", "coloc_with_susie", "results_all", str_c("coloc.susie_results_", results_names[i],".tsv")), sep = "\t", row.names = F, quote = F)
  
  }

}

