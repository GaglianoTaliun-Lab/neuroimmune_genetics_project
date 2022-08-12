# Description: run partial correlations for GWAS traits

# Load packages -----------------------------------------------------------

library(here)
library(LAVA)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

args <- commandArgs(TRUE)

# LAVA run version (v2, v2, v3):
run_version <- as.character(args[1])

# GWAS phenotypes to include in the partial correlation analysis:
gwas_phenos <- array(args[2:length(args)]) %>% sort()

args <-
  list(
    ref_prefix = here(project_dir,"reference_data","g1000_eur","g1000_eur"),
    loc_file = here(project_dir,"test_loci",stringr::str_c("gwas_filtered.loci.",run_version)),
    info_file = here(project_dir,"info_files","input.info.txt"),
    sample_overlap_file = here(project_dir,"sample_overlap",stringr::str_c("sample_overlap.",run_version,".txt")),
    phenotypes = gwas_phenos,
    output_filename = str_c(gwas_phenos, collapse = ":")
  )

print(args)

# name of LAVA bivar RDS file with results
bivar_file <- "AD_schwartzentruber2021:ALS_vanrheenen2021:CD_delange2017:LBD_chia2021:MS_imsgc2019:MSA_scholz2021:PD_nalls2019:SCZ_pardinas2018:UC_delange2017.bivar.lava.rds"

# Load data ---------------------------------------------------------------

loci <- LAVA::read.loci(args$loc_file)
input <-
  LAVA::process.input(
    input.info.file = args$info_file,
    sample.overlap.file = args$sample_overlap_file,
    ref.prefix = args$ref_prefix,
    phenos = args$phenotypes
  )

# Identify loci of interest -----------------------------------------------

bivar_results <-
  readRDS(here(project_dir, "lava_results","lava_GWAS",bivar_file)
  ) %>%
  purrr::discard(is.null) %>%
  qdapTools::list_df2df() %>%
  dplyr::select(-X1)

print(paste0("P-value threshold for bivar analysis: ", 0.05/nrow(bivar_results)))

# Filter for significant loci where:
# (i) > 2 distinct phenotypes
# (ii) all neuropsychiatric phenotypes are present
loci_of_interest <-
  bivar_results %>%
  dplyr::filter(p < 0.05/nrow(bivar_results)) %>%
  dplyr::group_by(locus) %>%
  dplyr::summarise(
    distinct_phen =
      list(c(phen1, phen2)) %>%
      lapply(., unique),
    n_distinct_phen =
      distinct_phen %>%
      unlist() %>%
      length(),
  ) %>%
  dplyr::filter(
    n_distinct_phen > 2
  )

print(paste0("Loci of interest = ", nrow(loci_of_interest)))

# Combinations
combn <-
  gwas_phenos %>%
  combn(m = 2) %>%
  t() %>%
  as_tibble() %>%
  dplyr::rename(
    phen1 = V1,
    phen2 = V2
  )

# Main --------------------------------------------------------------------

# Create empty lists for results
results <-
  vector(mode = "list", length = nrow(loci_of_interest)) %>%
  lapply(., function(x) x <- list())

for(i in 1:nrow(loci_of_interest)){

  locus_of_interest <-
    loci_of_interest %>%
    dplyr::slice(i)

  print(Sys.time())
  print(str_c("Locus: ", locus_of_interest$locus))

  # Process locus
  locus <-
    LAVA::process.locus(
      locus =
        loci %>%
        dplyr::filter(LOC == locus_of_interest$locus),
      input = input
    )

  # Extract some general locus info for the output
  loc_info <-
    data.frame(
      locus = locus$id,
      chr = locus$chr,
      start = locus$start,
      stop = locus$stop,
      n_snps = locus$n.snps,
      n_pcs = locus$K
    )

  for(j in 1:nrow(combn)){

    phen_to_corr <- c(combn$phen1[j], combn$phen2[j])
    phen_to_adj <- gwas_phenos[!gwas_phenos %in% phen_to_corr]

    print(Sys.time())
    print(str_c("Correlating: ", phen_to_corr[1], " ", phen_to_corr[2]))
    print(str_c("Adjusting for: ", phen_to_adj))

    # Partial corr
    loc_out <-
      run.pcor(locus, target = c(phen_to_corr), phenos = c(phen_to_adj))

    results[[i]][[j]] <-
      loc_info %>%
      dplyr::bind_cols(loc_out)

  }

}

# Save data ---------------------------------------------------------------
out_dir <- here(project_dir,"lava_results","lava_GWAS","partial_corr")

saveRDS(
  results,
  file = file.path(out_dir, str_c(args$output_filename, ".partcorr.lava.rds"))
)

print(str_c("Done! Analysis output written to ", out_dir, "/", args$output_filename, ".partcorr.lava.rds"))
