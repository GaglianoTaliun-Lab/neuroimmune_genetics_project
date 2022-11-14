library(dplyr)
library(here)
library(stringr)

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

###------------------------------------------------- arguments
args <- commandArgs(TRUE)

# prefix of LAVA output for that specific run:
datasets <- as.character(args[1])

# read RDS files
bivar <- readRDS(here(project_dir,"lava_results","lava_GWAS",stringr::str_c(datasets,".bivar.lava.rds")))
univar <- readRDS(here(project_dir,"lava_results","lava_GWAS",stringr::str_c(datasets,".univ.lava.rds")))

###---------------------------------------- Main

# Univariate -----------------------------------

# get number of univariate tests (number of loci in LAVA locus input)
ntest_univ <- read.table(here(project_dir, "test_loci", "gwas_filtered.loci_no_proxies"), sep = "\t", header = T) %>% 
  nrow(.)

# create dataframe with all univariate results
univ_all <- bind_rows(univar)

# create dataframe with only significant univariate results
univ_sign <- univ_all %>%
  filter(., p < 0.05/ntest_univ)

# Bivariate ------------------------------------

# get number of bivariate tests
ntests = 0
for (i in 1:length(bivar)) {
  condition = nrow(bivar[[i]]) > 0
  if (condition == "TRUE" && length(condition) != 0) {
    ntests <- ntests + nrow(bivar[[i]])
  }
}

# get pvalue bonferroni threshold
pvalue_bivar = 0.05/ntests

# remove null loci
bivar <- bivar[!sapply(bivar,is.null)]
bivar_all <- bind_rows(bivar)

# filter out non significant bivariate tests
bivar_significant <- lapply(bivar, function(x) filter(x, p <= pvalue_bivar))
bivar_significant <- bind_rows(bivar_significant)

# obtain number of significant bivariate tests
n_sig = nrow(bivar_significant)

### ---------------------------------------------- Save files

# write all results into table
write.table(bivar_all, here(project_dir,"lava_results","tsv_univ_bivar_GWAS","phenotypes_all_results_no_proxies_bivar.tsv"), sep = "\t", quote = F, row.names = F)
write.table(univ_all, here(project_dir,"lava_results","tsv_univ_bivar_GWAS","phenotypes_all_results_no_proxies_univ.tsv"), sep = "\t", quote = F, row.names = F)

# write significant results into table
write.table(bivar_significant, here(project_dir,"lava_results","tsv_univ_bivar_GWAS","phenotypes_significant_results_no_proxies_bivar.tsv"), sep = "\t", quote = F, row.names = F)
write.table(univ_sign, here(project_dir,"lava_results","tsv_univ_bivar_GWAS","phenotypes_significant_results_no_proxies_univ.tsv"), sep = "\t", quote = F, row.names = F)
