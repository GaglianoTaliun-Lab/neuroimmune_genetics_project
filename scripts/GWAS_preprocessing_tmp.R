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

# Load general data: gwas_filtered.loci from 100kb window:
loci <- read.table(paste0(project_dir,"/test_loci/window_100000/gene_filtered.loci"), sep = "\t", header = T)

out_dir <- paste0(project_dir,"/GWAS_summary_statistics/LAVA_sumstats/QTLs")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | FTD 2014 |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "GWAS_summary_statistics/FTD_GWAS_META.txt")
fread(file_path, nrows = 6)

FTD <- fread(file_path) %>%
  dplyr::mutate(CHR = as.factor(chr),
		N = 2154 + 4308) %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "a", "A") %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "g", "G") %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "c", "C") %>%
  dplyr::mutate_at(c("Allele1","Allele2"), stringr::str_replace, "t", "T") %>%
  dplyr::select(
    CHR,
    BP = Bp,
    A1 = Allele1,
    A2 = Allele2,
    BETA = beta1,
    SE,
    P = pValue,
    N) %>%
  dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
  colochelpR::convert_loc_to_rs(df = ., dbsnp_144)

# then need to remove positions with more than one rsid:
FTD <-
 FTD %>%
 dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
 dplyr::group_by(CHR_BP) %>%
 dplyr::filter(!any(row_number() > 1))  %>%
 dplyr::ungroup() %>%
 dplyr::select(-CHR_BP)

fwrite(FTD, file = here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","GWAS",stringr::str_c("FTD_ferrari2014", lava_ext)), sep = "\t")


