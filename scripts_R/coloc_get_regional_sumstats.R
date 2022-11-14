# create regional files for GWAS summary statistics to be used for coloc

library(here)
library(tidyverse)
library(stringr)
library(data.table)
library(rutils)
library(R.utils)
library(MafDb.1Kgenomes.phase3.hs37d5)
library(GenomicScores)
library(colochelpR)

# Arguments ------------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

mafdb <- MafDb.1Kgenomes.phase3.hs37d5

GWAS_coloc = c("AD_schwartzentruber2021.lava.gz", "CD_delange2017.lava.gz", "ALS_vanrheenen2021.lava.gz", "MS_imsgc2019.lava.gz",
               "PD_nalls2019.lava.gz", "SCZ_pardinas2018.lava.gz", "UC_delange2017.lava.gz")

mafs <- fread(here(project_dir, "reference_data", "g1000_eur", "g1000_eur.frq")) %>% dplyr::select(SNP, MAF)

# Load files -----------------------------------------------------------

gene_table <- read.table(here(project_dir, "colocalization", "gene_table_coloc_FDR.txt"),
                        sep = "\t", header = T)

lava_sumstats_files <- list.files(here(project_dir, "GWAS_summary_statistics", "LAVA_sumstats", "GWAS"), 
                                  pattern = ".lava.gz", full.names = T)

# Format dataframes and save regional sumstats -----------------------------------------------------------

for (i in 1:length(lava_sumstats_files)) {
  
  if (basename(lava_sumstats_files[i]) %in% GWAS_coloc) {
    
    curr_trait = basename(lava_sumstats_files[i]) %>% stringr::str_remove(., "_[:graph:]+")
    gene_table_trait <- filter(gene_table, trait == curr_trait)
    
    lava_sumstats <- fread(lava_sumstats_files[i])
    
    for (j in 1:nrow(gene_table_trait)) {
      
      chr_pos = gene_table_trait$chr[j]
      start_pos = gene_table_trait$start[j]
      end_pos = gene_table_trait$end[j]
      
      if ("BETA" %in% names(lava_sumstats)) {
        
        lava_sumstats_region <- filter(lava_sumstats, CHR == chr_pos & BP >= start_pos & BP <= end_pos) %>%
          mutate(GWAS = curr_trait) %>%
          dplyr::select(., GWAS, SNP, CHR, BP, beta = BETA, se = SE, pvalues = P, A1, A2, N) %>%
          colochelpR::get_varbeta(.)
        
      } else if ("OR" %in% names(lava_sumstats)) {
        
        lava_sumstats_region <- filter(lava_sumstats, CHR == chr_pos & BP >= start_pos & BP <= end_pos) %>%
          mutate(GWAS = curr_trait) %>%
          dplyr::select(., GWAS, SNP, CHR, BP, OR, pvalues = P, A1, A2, N)

      } else next
      
#      mafs <- GenomicScores::gscores(x = mafdb, ranges = unique(lava_sumstats_region$SNP) %>% as.character(), pop = "EUR_AF")
#      mafs <- mafs %>%
#        as.data.frame() %>%
#        tibble::rownames_to_column(var = "SNP") %>%
#        dplyr::rename(MAF = EUR_AF) %>%
#        dplyr::select(SNP, MAF)
      
#      head(mafs, 5)
      
      lava_sumstats_region <- lava_sumstats_region %>%
        inner_join(mafs, by = "SNP")
      
      head(lava_sumstats_region, 5)
      
      write.table(lava_sumstats_region, here(project_dir, "colocalization", "GWAS_regional_sumstats",
                                             stringr::str_c(curr_trait, "_", gene_table_trait$gene[j], ".tsv")),
                  sep = "\t", row.names = F, quote = F)
      
    }
    
  } else next
  
}
