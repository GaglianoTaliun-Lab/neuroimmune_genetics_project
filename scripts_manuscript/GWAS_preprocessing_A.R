# process GWAS summary statistics for LAVA - to run as a normal job on compute canada

# Summary statistics must have the following columns:
# SNP / ID / SNPID_UKB/ SNPID / MarkerName / RSID / RSID_UKB: SNP IDs
# A1 / ALT: effect allele
# A2 / REF: reference allele
# N / NMISS / OBS_CT / N_analyzed: number of samples
# Z / T / STAT / Zscore: if provided, no p-values or coefficients are needed; otherwise, please provide both:
# B / BETA / OR / logOdds: effect size coefficients
# P: p-values

# install.packages("devtools")
# install.packages("BiocManager")

# BiocManager::install("biomaRt")
# BiocManager::install("BSgenome")

library(here)
library(devtools)
library(BiocManager)
library(dplyr)
library(data.table)
library(biomaRt)
library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

# devtools::install_github("RHReynolds/colochelpR")
# devtools::install_github("RHReynolds/rutils")
library(colochelpR)
library(rutils)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"
lava_ext = ".lava.gz"
dbsnp_144 <- SNPlocs.Hsapiens.dbSNP144.GRCh37
path_hg38_hg37 <- here(project_dir,"reference_data","hg38ToHg19.over.chain")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | PD |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path <- file.path(project_dir, "GWAS_summary_statistics/PD_nalls2019.tab")
fread(file_path, nrows = 6)

PD <- fread(file_path) %>%
 filter(., freq >= 0.01) %>%
 dplyr::mutate(CHR = as.factor(chr)) %>%
 dplyr::mutate(N = N_cases + N_controls) %>%
 dplyr::mutate(N = N + 0.0) %>%
 dplyr::select(
   CHR,
   BP = pos,
   A1,
   A2,
   BETA = b,
   SE = se,
   P = p,
   N) %>%
  dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
  colochelpR::convert_loc_to_rs(df = ., dbsnp_144)

# then need to remove positions with more than one rsid:
PD <-
 PD %>%
 dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
 dplyr::group_by(CHR_BP) %>%
 dplyr::filter(!any(row_number() > 1))  %>%
 dplyr::ungroup() %>%
 dplyr::select(-CHR_BP)

fwrite(PD, file = here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","GWAS",stringr::str_c("PD_nalls2019", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | AD_wightman |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path <- file.path(project_dir, "GWAS_summary_statistics/AD_wightman2021.txt")
fread(file_path, nrows = 6)

AD_w <- fread(file_path) %>%
 dplyr::mutate(CHR = as.factor(chr)) %>%
 dplyr::mutate(N = N + 0.0) %>%
 dplyr::select(
   CHR,
   BP = PosGRCh37,
   A1 = testedAllele,
   A2 = otherAllele,
   Z = z,
   P = p,
   N) %>%
  dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
  colochelpR::convert_loc_to_rs(df = ., dbsnp_144)

# then need to remove positions with more than one rsid:
AD_w <-
 AD_w %>%
 dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
 dplyr::group_by(CHR_BP) %>%
 dplyr::filter(!any(row_number() > 1))  %>%
 dplyr::ungroup() %>%
 dplyr::select(-CHR_BP)

fwrite(AD_w, file = here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","GWAS",stringr::str_c("AD_wightman2021", lava_ext)), sep = "\t")
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | AD_schwartzentruber |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "GWAS_summary_statistics/AD_schwartzentruber2021.tsv")
fread(file_path, nrows = 6)

AD_s <- fread(file_path) %>%
 dplyr::mutate(N = 75671 + 397844, CHR = as.factor(chromosome)) %>%
 dplyr::select(
   SNP = variant_id,
   CHR,
   BP = base_pair_location,
   A1 = effect_allele,
   A2 = other_allele,
   BETA = beta,
   SE = standard_error,
   P = p_value,
   N) %>%
  dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":"))

fwrite(AD_s, file = here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","GWAS",stringr::str_c("AD_schwartzentruber2021", lava_ext), sep = "\t")
