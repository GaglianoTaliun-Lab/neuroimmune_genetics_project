library(tidyverse)
library(stringr)
library(data.table)
library(here)

project_dir="/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

pd_raw <- read.table(here(project_dir,"GWAS_summary_statistics","GWAS_compressed_summary_statistics","PD_data","PD_blauwendraat2019.tbl"), sep = "\t", header = T)

# pd_test <- pd_raw[c(1:1000),]

# Obtain CHR and BP columns, put alleles as big letters:

pd_filtered <- pd_raw %>%
  tidyr::separate(., col = MarkerName, into = c("CHR","BP")) %>%
  dplyr::mutate_at("CHR", stringr::str_replace, "chr", "") %>%
  dplyr::mutate(CHR = as.factor(CHR)) %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "a", "A") %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "g", "G") %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "c", "C") %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "t", "T") %>%
# Filtering according to Cornelis' suggestions:
  filter(., Freq1 > 0.01) %>%
  filter(., HetISq < 80) %>%
  filter(., HetDf >= 8)

# Obtain N for each variant:
samples <- as.data.frame(pd_filtered$Direction)
colnames(samples) <- "Direction"
samples <- samples %>% 
  tidyr::separate(., Direction, into = c("s0","s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12","s13"), sep = "") %>%
  dplyr::mutate(N = 0) %>% .[,c(2:15)]

N_per_cohort <- read.table(here(project_dir,"GWAS_summary_statistics","GWAS_compressed_summary_statistics","PD_data","sample_size_per_cohort.txt"), sep="\t", header = T) %>%
  dplyr::mutate(Cases = as.numeric(Cases), Controls = as.numeric(Controls), Total = as.numeric(Total))

samples$N <- mapply(function(x, y) if (x != "?") y = y + N_per_cohort[1,4] else y = y, samples$s1, samples$N)
samples$N <- mapply(function(x, y) if (x != "?") y = y + N_per_cohort[2,4] else y = y, samples$s2, samples$N)
samples$N <- mapply(function(x, y) if (x != "?") y = y + N_per_cohort[3,4] else y = y, samples$s3, samples$N)
samples$N <- mapply(function(x, y) if (x != "?") y = y + N_per_cohort[4,4] else y = y, samples$s4, samples$N)
samples$N <- mapply(function(x, y) if (x != "?") y = y + N_per_cohort[5,4] else y = y, samples$s5, samples$N)
samples$N <- mapply(function(x, y) if (x != "?") y = y + N_per_cohort[6,4] else y = y, samples$s6, samples$N)
samples$N <- mapply(function(x, y) if (x != "?") y = y + N_per_cohort[7,4] else y = y, samples$s7, samples$N)
samples$N <- mapply(function(x, y) if (x != "?") y = y + N_per_cohort[8,4] else y = y, samples$s8, samples$N)
samples$N <- mapply(function(x, y) if (x != "?") y = y + N_per_cohort[9,4] else y = y, samples$s9, samples$N)
samples$N <- mapply(function(x, y) if (x != "?") y = y + N_per_cohort[10,4] else y = y, samples$s10, samples$N)
samples$N <- mapply(function(x, y) if (x != "?") y = y + N_per_cohort[11,4] else y = y, samples$s11, samples$N)
samples$N <- mapply(function(x, y) if (x != "?") y = y + N_per_cohort[12,4] else y = y, samples$s12, samples$N)
samples$N <- mapply(function(x, y) if (x != "?") y = y + N_per_cohort[13,4] else y = y, samples$s13, samples$N)

pd_final <- cbind(pd_filtered, N = samples$N)
write.table(pd_final, here(project_dir,"GWAS_summary_statistics","PD_blauwendraat2019_filtered.tbl", sep = "\t", row.names = F, quote = F)

