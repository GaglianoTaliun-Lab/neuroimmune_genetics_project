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
#                                               | CD |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "GWAS_summary_statistics/CD_delange2017.tsv")
fread(file_path, nrows = 6)

CD <- fread(file_path) %>%
  dplyr::mutate(N = 12194 + 28072, CHR = as.factor(chromosome)) %>%
  dplyr::select(
    CHR,
    BP = base_pair_location,
    A1 = effect_allele, 
    A2 = other_allele, 
    BETA = beta,
    SE = standard_error,
    P = p_value,
    N) %>%
  dplyr::mutate(A1 = stringr::str_to_upper(A1), A2 = stringr::str_to_upper(A2)) %>%
  dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
  colochelpR::convert_loc_to_rs(df = ., dbsnp_144)

# then need to remove positions with more than one rsid:
CD <-
  CD %>%
  dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
  dplyr::group_by(CHR_BP) %>%
  dplyr::filter(!any(row_number() > 1))  %>%
  dplyr::ungroup() %>%
  dplyr::select(-CHR_BP)

fwrite(CD, file = here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","GWAS",stringr::str_c("CD_delange2017", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | UC |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "GWAS_summary_statistics/UC_delange2017.tsv")
fread(file_path, nrows = 6)

UC <- fread(file_path) %>%
  dplyr::mutate(N = 12366 + 33609, CHR = as.factor(chromosome)) %>%
  dplyr::select(
    CHR,
    BP = base_pair_location,
    A1 = effect_allele, 
    A2 = other_allele, 
    BETA = beta,
    SE = standard_error,
    P = p_value,
    N) %>%
  dplyr::mutate(A1 = stringr::str_to_upper(A1), A2 = stringr::str_to_upper(A2)) %>%
  dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
  colochelpR::convert_loc_to_rs(df = ., dbsnp_144)

# then need to remove positions with more than one rsid:
UC <-
  UC %>%
  dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
  dplyr::group_by(CHR_BP) %>%
  dplyr::filter(!any(row_number() > 1))  %>%
  dplyr::ungroup() %>%
  dplyr::select(-CHR_BP)

fwrite(UC, file = here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","GWAS",stringr::str_c("UC_delange2017", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | MS |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "GWAS_summary_statistics/MS_sawcer2011.tsv")
fread(file_path, nrows = 6)

MS <- fread(file_path) %>%
  dplyr::mutate(N = 9772.0 + 17376.0, CHR = as.factor(chromosome)) %>%
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

MS2 <- MS[!grepl("X", MS$CHR),]

fwrite(MS2, file = here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","GWAS",stringr::str_c("MS_sawcer2011", lava_ext)), sep = "\t")

head(MS)
