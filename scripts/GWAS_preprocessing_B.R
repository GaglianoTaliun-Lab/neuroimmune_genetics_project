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
#                                               | LBD |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "GWAS_summary_statistics/LBD_chia2021.tsv")
fread(file_path, nrows = 6)

LBD <- fread(file_path) %>%
 dplyr::mutate(N = 2591 + 4027, CHR = as.factor(chromosome)) %>%
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
  dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
 rutils::liftover_coord(df = ., path_to_chain = path_hg38_hg37)

fwrite(LBD, file = here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","GWAS",stringr::str_c("LBD_chia2021", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | ALS 2016 |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "GWAS_summary_statistics/ALS_vanrheenen2016.txt")
fread(file_path, nrows = 6)

ALS1 <- fread(file_path) %>%
  dplyr::mutate(N = 12577 + 23475, CHR = as.factor(chr)) %>%
  dplyr::select(
    SNP = snp,
    CHR,
    BP = bp,
    A1 = a1,
    A2 = a2,
    BETA = b,
    SE = se,
    P = p,
    N) %>%
  dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":"))

fwrite(ALS1, file = here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","GWAS",stringr::str_c("ALS_vanrheenen2016", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | SCZ |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "GWAS_summary_statistics/SCZ_pardinas2018.txt")
fread(file_path, nrows = 6)

SCZ <- fread(file_path) %>%
  dplyr::mutate(N = 105318, CHR = as.factor(CHR), BETA = log(OR)) %>%
  dplyr::select(
    SNP,
    CHR,
    BP,
    A1,
    A2,
    BETA,
    SE,
    P,
    N) %>%
  dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":"))

fwrite(SCZ, file = here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","GWAS",stringr::str_c("SCZ_pardinas2018", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | ALS 2021 |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "GWAS_summary_statistics/ALS_vanrheenen2021.txt")
fread(file_path, nrows = 6)

ALS2 <- fread(file_path) %>%
  dplyr::mutate(CHR = as.factor(chromosome)) %>%
  dplyr::mutate_at(c("effect_allele","other_allele"), stringr::str_replace, "a", "A") %>%
  dplyr::mutate_at(c("effect_allele","other_allele"), stringr::str_replace, "g", "G") %>%
  dplyr::mutate_at(c("effect_allele","other_allele"), stringr::str_replace, "c", "C") %>%
  dplyr::mutate_at(c("effect_allele","other_allele"), stringr::str_replace, "t", "T") %>%
  dplyr::select(
    SNP = rsid,
    CHR,
    BP = base_pair_location,
    A1 = effect_allele,
    A2 = other_allele,
    BETA = beta,
    SE = standard_error,
    P = p_value,
    N = N_effective) %>%
  dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":"))

fwrite(ALS2, file = here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","GWAS",stringr::str_c("ALS_vanrheenen2021", lava_ext)), sep = "\t")


