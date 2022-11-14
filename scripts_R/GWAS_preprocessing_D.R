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
library(stringr)
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
#                                               | PD_blauwendraat |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path <- file.path(project_dir, "GWAS_summary_statistics/PD_blauwendraat2019_filtered.tbl")
fread(file_path, nrows = 6)

PD <- fread(file_path) %>%
 dplyr::mutate(CHR = as.factor(CHR)) %>%
 dplyr::select(
   CHR,
   BP,
   A1 = Allele1,
   A2 = Allele2,
   BETA = Effect,
   SE = StdErr,
   P = P.value,
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

fwrite(PD, file = here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","GWAS",stringr::str_c("PD_blauwendraat2019", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | MS_imsgc |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path <- file.path(project_dir, "GWAS_summary_statistics/MS_imsgc2019.meta")
fread(file_path, nrows = 6)

MS <- fread(file_path) %>%
 dplyr::mutate(CHR = as.factor(CHR), N = 14802 + 26703) %>%
 dplyr::select(
   CHR,
   BP,
   A1,
   A2,
   OR = OR,
   P = P,
   N) %>%
  dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
  colochelpR::convert_loc_to_rs(df = ., dbsnp_144)

# then need to remove positions with more than one rsid:
MS <-
 MS %>%
 dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
 dplyr::group_by(CHR_BP) %>%
 dplyr::filter(!any(row_number() > 1))  %>%
 dplyr::ungroup() %>%
 dplyr::select(-CHR_BP)

fwrite(MS, file = here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","GWAS",stringr::str_c("MS_imsgc2019", lava_ext)), sep = "\t")
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | AD_kunkle |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "GWAS_summary_statistics/AD_kunkle2019.txt")
fread(file_path, nrows = 6)

AD_k <- fread(file_path) %>%
 dplyr::mutate(N = 21982 + 41944, CHR = as.factor(Chromosome)) %>%
 dplyr::select(
   SNP = MarkerName,
   CHR,
   BP = Position,
   A1 = Effect_allele,
   A2 = Non_Effect_allele,
   BETA = Beta,
   SE,
   P = Pvalue,
   N) %>%
  dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":"))

fwrite(AD_k, file = here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","GWAS",stringr::str_c("AD_kunkle2019", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | MSA |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path <- file.path(project_dir, "GWAS_summary_statistics/MSA_scholz2021.txt")
fread(file_path, nrows = 6)

MSA <- fread(file_path) %>%
 dplyr::mutate(CHR = as.factor(CHROM), N = CASE_ALLELE_CT + CTRL_ALLELE_CT) %>%
 dplyr::select(
   CHR,
   BP = POS,
   A1,
   A2,
   BETA,
   SE,
   P,
   N) %>%
  dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
 rutils::liftover_coord(df = ., path_to_chain = path_hg38_hg37) %>%
  colochelpR::convert_loc_to_rs(df = ., dbsnp_144)

# then need to remove positions with more than one rsid:
MSA <-
 MSA %>%
 dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
 dplyr::group_by(CHR_BP) %>%
 dplyr::filter(!any(row_number() > 1))  %>%
 dplyr::ungroup() %>%
 dplyr::select(-CHR_BP)

fwrite(MSA, file = here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","GWAS",stringr::str_c("MSA_scholz2021", lava_ext)), sep = "\t")
