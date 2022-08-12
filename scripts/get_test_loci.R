# Description: get test loci from GWAS traits

# Load packages -----------------------------------------------------------

library(here)
library(data.table)
library(GenomicRanges)
library(qdapTools)
library(stringr)
library(tidyverse)
library(rutils)
library(R.utils)

# Arguments ---------------------------------------------------------------
args <- commandArgs(TRUE)

# LAVA run version (v1, v2, v3):
run_version <- as.character(args[1])

# Load data ---------------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# LD blocks obtained from: https://github.com/cadeleeuw/lava-partitioning/raw/main/LAVA_s2500_m25_f1_w200.blocks
ld_blocks <-
  read.table(here(project_dir,"reference_data","blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"),
             sep = " ", header = T)

gwas_list <-
  setNames(
    object = list(
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","AD_schwartzentruber2021.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","ALS_vanrheenen2021.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","CD_delange2017.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","LBD_chia2021.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","MS_imsgc2019.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","MSA_scholz2021.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","PD_nalls2019.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","SCZ_pardinas2018.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","UC_delange2017.lava.gz"))
    ),
    nm = c("AD_schwartzentruber2021", "ALS_vanrheenen2021", "CD_delange2017", "LBD_chia2021", "MS_imsgc2019", "MSA_scholz2021", "PD_nalls2019", "SCZ_pardinas2018", "UC_delange2017")
  )

# Main --------------------------------------------------------------------

############### !!! fread reads pvalues for Kunkle et al. 2019 as chr, leading to bad filtering of pvalues below.
# need to investigate why, but this is a quick fix:

# gwas_list$AD_kunkle2019[[8]] <- as.numeric(gwas_list$AD_kunkle2019[[8]])

# Filter for only genome-wide significant loci and convert to granges
gr_list <-
  gwas_list %>%
  lapply(., function(gwas){
    gwas %>%
      dplyr::filter(P < 5e-8) %>%
      GenomicRanges::makeGRangesFromDataFrame(
        .,
        keep.extra.columns = TRUE,
        ignore.strand = TRUE,
        seqinfo = NULL,
        seqnames.field = "CHR",
        start.field = "BP",
        end.field = "BP"
      )
  })

# Add locus id and rename remaining columns to fit LAVA requirements
ld_blocks <-
  ld_blocks %>%
  dplyr::rename_with(
    .fn = stringr::str_to_upper,
    .cols = everything()
  ) %>%
  dplyr::mutate(
    LOC = dplyr::row_number()
  ) %>%
  dplyr::select(LOC, everything())

# Convert to granges
ld_blocks_gr <-
  ld_blocks %>%
  GenomicRanges::makeGRangesFromDataFrame(
    .,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqinfo = NULL,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "stop"
  )

# Overlap granges objects
overlap_list <-
  gr_list %>%
  lapply(., function(gr){
    GenomicRanges::findOverlaps(gr, ld_blocks_gr, type = "within") %>%
      as_tibble()
  })


# Extract relevant rows from ld_blocks using overlap indices
# Rename remaining columns to fit with LAVA requirements
loci <-
  ld_blocks %>%
  dplyr::slice(
    overlap_list %>%
      qdapTools::list_df2df(col1 = "gwas") %>%
      .[["subjectHits"]] %>%
      unique()
    ) %>%
  dplyr::arrange(LOC)

# Generate df of GWAS loci that overlap which LD blocks
overlap_df_list <- vector(mode = "list", length = length(gr_list))

for(i in 1:length(gr_list)){

  gr <- gr_list[[i]]

  names(overlap_df_list)[i] <- names(gr_list)[i]

  if(names(gr_list)[i] == c("MS_imsgc2019")){
    
    overlap_df_list[[i]] <-
      gr %>%
      as_tibble() %>%
      dplyr::mutate(
        BETA = NA,
        SE = NA
      ) %>%
      dplyr::select(seqnames, start, end, SNP, A1, A2, BETA, SE, OR, P, N)

  } else{
    
    overlap_df_list[[i]] <-
      gr %>%
      as_tibble() %>%
      dplyr::mutate(
        OR = NA,
      ) %>%
      dplyr::select(seqnames, start, end, SNP, A1, A2, BETA, SE, OR, P, N)
  }
  
  overlap_df_list[[i]] <-
    overlap_df_list[[i]] %>%
    dplyr::rename_with(
      ~ stringr::str_c("GWAS", .x, sep = "_")
    ) %>%
    dplyr::slice(overlap_list[[i]]$queryHits) %>%
    dplyr::bind_cols(
      ld_blocks_gr %>%
        as_tibble() %>%
        dplyr::rename_with(
          ~ stringr::str_c("LD", .x, sep = "_")
        ) %>%
        dplyr::slice(overlap_list[[i]]$subjectHits)
    ) %>%
    dplyr::select(-contains("strand")) %>%
    dplyr::rename_with(
      ~ stringr::str_replace(.x,
                             pattern = "seqnames",
                             replacement = "CHR"),
      .col = dplyr::ends_with("seqnames"))
}

overlap_df <-
  overlap_df_list %>%
  qdapTools::list_df2df(col1= "GWAS")

# Save data ---------------------------------------------------------------

write.table(
  loci,
  here(project_dir,"test_loci", stringr::str_c("gwas_filtered.loci.", run_version)),
  sep = "\t",
  row.names = F,
  quote = F
  )

saveRDS(
  overlap_df,
  file = here(project_dir,"test_loci", stringr::str_c("gwas_filtered_loci." ,run_version, ".rds"))
  )
