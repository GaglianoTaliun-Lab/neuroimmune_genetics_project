# check the correlation between MAFs from original datasets and MAFs from reference (1KGP-EUR)

library(here)
library(tidyverse)
library(stringr)
library(data.table)
library(rutils)
library(R.utils)
library(MafDb.1Kgenomes.phase3.hs37d5)
library(GenomicScores)

# Arguments ------------------------------------------------------------
project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"
GWAS_coloc_traits = c("AD", "ALS", "PD")

coloc_files <- list()
coloc_dfs <- list()

# Load and format files -----------------------------------------------------------

# load coloc regional summary statistics:
for (i in 1:length(GWAS_coloc_traits)) {
  
  coloc_files[[i]] <- list.files(path = here(project_dir, "colocalization", "GWAS_regional_sumstats"), 
                                 pattern = stringr::str_c(GWAS_coloc_traits[i],"_*"), full.names = T) %>%
    lapply(., function(x) {fread(x)})
  
  coloc_dfs[[i]] <- rbindlist(coloc_files[[i]])
  
}

names(coloc_dfs) <- GWAS_coloc_traits

# merging MAFs from original datasets

GWAS_sumstats <- list()
# ----------------------------------------------------------- AD
GWAS_sumstats[[1]] <- read.table(here(project_dir, "GWAS_summary_statistics", "AD_schwartzentruber2021.tsv"), sep = "\t", header = T) %>%
  dplyr::select(SNP = variant_id,
                AF = effect_allele_frequency) %>%
  mutate(MAF_GWAS = case_when(
    AF > 0.5 ~ 1 - AF,
    AF < 0.5 ~ AF,
    AF == 0.5 ~ AF
  )) %>%
  right_join(., coloc_dfs$AD, by = "SNP")
  
AD_cor <- cor.test(GWAS_sumstats[[1]]$MAF, GWAS_sumstats[[1]]$MAF_GWAS)

# ----------------------------------------------------------- ALS
GWAS_sumstats[[2]] <- read.table(here(project_dir, "GWAS_summary_statistics", "ALS_vanrheenen2021.txt"), sep = "\t", header = T) %>%
  mutate(AF = as.numeric(effect_allele_frequency)) %>%
  dplyr::select(SNP = rsid,
                AF) %>%
  mutate(MAF_GWAS = case_when(
    AF > 0.5 ~ 1 - AF,
    AF < 0.5 ~ AF,
    AF == 0.5 ~ AF
  )) %>%
  right_join(., coloc_dfs$ALS, by = "SNP")

ALS_cor <- cor.test(GWAS_sumstats[[2]]$MAF, GWAS_sumstats[[2]]$MAF_GWAS)

# ----------------------------------------------------------- PD
GWAS_sumstats[[3]] <- read.table(here(project_dir, "GWAS_summary_statistics", "PD_nalls2019.tab"), sep = "\t", header = T) %>%
  mutate(AF = as.numeric(freq)) %>%
  dplyr::select(CHR = chr,
                BP = pos,
                AF) %>%
  mutate(MAF_GWAS = case_when(
    AF > 0.5 ~ 1 - AF,
    AF < 0.5 ~ AF,
    AF == 0.5 ~ AF
  )) %>%
  right_join(., coloc_dfs$PD, by = c("CHR", "BP"))

PD_cor <- cor.test(GWAS_sumstats[[3]]$MAF, GWAS_sumstats[[3]]$MAF_GWAS)

# ----------------------------------------------------------- 

# merge dfs into one large df:
GWAS_sumstats_df <- rbindlist(GWAS_sumstats, use.names = TRUE)
global <- cor.test(GWAS_sumstats_df$MAF, GWAS_sumstats_df$MAF_GWAS)
global

# save all cor.test results in table
data.frame("trait" = c(GWAS_coloc_traits, "global"),
           "cor_estimate" = c(AD_cor$estimate, ALS_cor$estimate, PD_cor$estimate, global$estimate),
           "p_value" = c(AD_cor$p.value, ALS_cor$p.value, PD_cor$p.value, global$p.value)) %>%
  write.table(., here(project_dir, "colocalization", "corr_results_MAF_comparison.tsv"), sep = "\t", row.names = F, quote = F)

# plot in scatterplot and save:
ggplot(GWAS_sumstats_df, aes(x = MAF, y = MAF_GWAS, color = GWAS)) +
  geom_point() +
  theme_bw()
ggsave(here(project_dir, "colocalization", "corr_results_MAF_comparison.pdf"))

