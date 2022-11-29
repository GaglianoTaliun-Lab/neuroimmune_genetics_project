# Description: script to obtain heat map plots of the global genetic correlations across traits.

# Libraries ----------------------------------------------------------------------------------------

library(tidyverse)
library(here)
library(gghighlight)
library(stringr)
library(reshape2)
library(corrplot)

# Arguments ----------------------------------------------------------------------------------------

project_dir <- "Documents/research-projects/neuro_immune/ldsc_results"

correlations <- read.table(here(project_dir, "ldsc_correlations.txt"), sep = "\t", header = T)

p_nominal = 0.05

source(here(project_dir, "global_corr_function_RHR.R"))

# Wrangle data ----------------------------------------------------------------------------------------

# round rg to 2 digits
correlations <- correlations %>%
  mutate(rg = round(rg, digits = 2))

# remove unwanted phenotypes and remove authors/year:
correlations_prox <- correlations %>% 
  filter(p1 != "AD_wightman2021") %>%
  filter(p2 != "AD_wightman2021") %>%
  filter(p1 != "AD_kunkle2019") %>%
  filter(p2 != "AD_kunkle2019") %>%
  filter(p1 != "SCZ_lam2019") %>%
  filter(p2 != "SCZ_lam2019") %>%
  filter(p1 != "MS_sawcer2011") %>%
  filter(p2 != "MS_sawcer2011") %>%
  filter(p1 != "PD_blauwendraat2019") %>%
  filter(p2 != "PD_blauwendraat2019") %>%
  filter(p1 != "ALS_vanrheenen2016") %>%
  filter(p2 != "ALS_vanrheenen2016") %>%
  filter(p1 != "MSA_scholz2021") %>%
  filter(p2 != "MSA_scholz2021") %>%
  mutate(p1 = str_remove(p1, "_.*"), p2 = str_remove(p2, "_.*"))

# remove unwanted phenotypes and remove authors/year:
correlations_noprox <- correlations %>% 
  filter(p1 != "AD_wightman2021") %>%
  filter(p2 != "AD_wightman2021") %>%
  filter(p1 != "AD_schwartzentruber2021") %>%
  filter(p2 != "AD_schwartzentruber2021") %>%
  filter(p1 != "SCZ_lam2019") %>%
  filter(p2 != "SCZ_lam2019") %>%
  filter(p1 != "MS_sawcer2011") %>%
  filter(p2 != "MS_sawcer2011") %>%
  filter(p1 != "PD_nalls2019") %>%
  filter(p2 != "PD_nalls2019") %>%
  filter(p1 != "ALS_vanrheenen2016") %>%
  filter(p2 != "ALS_vanrheenen2016") %>%
  filter(p1 != "MSA_scholz2021") %>%
  filter(p2 != "MSA_scholz2021") %>%
  filter(p1 != "NA") %>%
  mutate(p1 = str_remove(p1, "_.*"), p2 = str_remove(p2, "_.*"))

correlations_list <- list(correlations_prox, correlations_noprox)

# Using RHR function --------------------------------------------------------------------------------------------

correlations_prox2 <- correlations_prox %>%
  mutate(ord_p1 = factor(p1, levels = c("AD","LBD","PD","FTD","ALS","SCZ","MS","UC","CD")),
         ord_p2 = factor(p2, levels = c("AD","LBD","PD","FTD","ALS","SCZ","MS","UC","CD"))) %>%
  arrange(ord_p1, as.factor(ord_p1), as.factor(ord_p2))

plot_global_corr(global_corr = correlations_prox2, n_phenotypes = 9)

ggsave(here(project_dir,"global_correlations_significant_with_proxies_half_matrix.jpg"), width = 40, height = 30, units = "cm")

# --- without proxies:
correlations_prox2 <- correlations_noprox %>%
  mutate(ord_p1 = factor(p1, levels = c("AD","LBD","PD","FTD","ALS","SCZ","MS","UC","CD")),
         ord_p2 = factor(p2, levels = c("AD","LBD","PD","FTD","ALS","SCZ","MS","UC","CD"))) %>%
  arrange(ord_p1, as.factor(ord_p1), as.factor(ord_p2))

plot_global_corr(global_corr = correlations_prox2, n_phenotypes = 9)

ggsave(here(project_dir,"global_correlations_significant_no_proxies_half_matrix.jpg"), width = 40, height = 30, units = "cm")

# export genome-wide correlations table
correlations_prox %>%
  dplyr::select(Trait_1 = p1,
                Trait_2 = p2,
                rg,
                SE = se,
                Z_score = z,
                pvalue = p) %>%
write.table(., here(project_dir, "genome_wide_genetic_correlations_table.txt"), quote = F, sep = "\t", row.names = F)

correlations_noprox %>%
  dplyr::select(Trait_1 = p1,
                Trait_2 = p2,
                rg,
                SE = se,
                Z_score = z,
                pvalue = p) %>%
  write.table(., here(project_dir, "genome_wide_genetic_correlations_table_noproxies.txt"), quote = F, sep = "\t", row.names = F)

