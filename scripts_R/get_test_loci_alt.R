# Description: get test loci from GWAS traits
# script to automate the inclusion of Z, BETA, OR, SE columns for GWAS sumstats overlap (in comparison to 'get_test_loci.R'),
# also included the ALS_vanrheenen2021 sumstats instead of the 2016.

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

# Load data ---------------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# LD blocks obtained from: https://github.com/cadeleeuw/lava-partitioning/raw/main/LAVA_s2500_m25_f1_w200.blocks
ld_blocks <-
  read.table(here(project_dir,"reference_data","blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"),
             sep = " ", header = T)

# variable gwas list depending on the gwas sumstats to include (with proxies):
gwas_list_prox <-
  setNames(
    object = list(
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","AD_schwartzentruber2021.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","ALS_vanrheenen2021.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","CD_delange2017.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","LBD_chia2021.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","MS_imsgc2019.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","FTD_ferrari2014.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","PD_nalls2019.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","SCZ_pardinas2018.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","UC_delange2017.lava.gz"))
    ),
    nm = c("AD_schwartzentruber2021", "ALS_vanrheenen2021", "CD_delange2017", "LBD_chia2021", "MS_imsgc2019", "FTD_ferrari2014", "PD_nalls2019", "SCZ_pardinas2018", "UC_delange2017")
  )

# variable gwas list (without proxies):
gwas_list_noprox <-
  setNames(
    object = list(
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","AD_kunkle2019.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","ALS_vanrheenen2021.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","CD_delange2017.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","LBD_chia2021.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","MS_imsgc2019.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","FTD_ferrari2014.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","PD_blauwendraat2019.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","SCZ_pardinas2018.lava.gz")),
      fread(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","GWAS","UC_delange2017.lava.gz"))
    ),
    nm = c("AD_kunkle2019", "ALS_vanrheenen2021", "CD_delange2017", "LBD_chia2021", "MS_imsgc2019", "FTD_ferrari2014", "PD_blauwendraat2019", "SCZ_pardinas2018", "UC_delange2017")
  )

# Main --------------------------------------------------------------------

list_prox_noprox <- list(gwas_list_prox, gwas_list_noprox)

for (k in 1:length(list_prox_noprox)){
  
  # Filter for only genome-wide significant loci and convert to granges
  gr_list <-
    list_prox_noprox[[k]] %>%
    lapply(., function(each_gwas){
      each_gwas %>%
	mutate(P = as.numeric(P)) %>%
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
  # here I need to standardize the column names across gwas sumstats,
  # so that all sumstats have the same column names (e.g. add beta or Z columns filled with NA)
  overlap_df_list <- vector(mode = "list", length = length(gr_list))
  
  for(i in 1:length(gr_list)){
    
    gr <- gr_list[[i]]
    
    names(overlap_df_list)[i] <- names(gr_list)[i]
    
    if("OR" %in% names(list_prox_noprox[[k]][[i]])){
      
      overlap_df_list[[i]] <-
        gr %>%
        as_tibble() %>%
        dplyr::mutate(
          BETA = NA,
          SE = NA,
          Z = NA,
        ) %>%
        dplyr::select(seqnames, start, end, SNP, A1, A2, BETA, SE, OR, Z, P, N)
      
    } else if ("Z" %in% names(list_prox_noprox[[k]][[i]])){
      
      overlap_df_list[[i]] <-
        gr %>%
        as_tibble() %>%
        dplyr::mutate(
          BETA = NA,
          SE = NA,
          OR = NA,
        ) %>%
        dplyr::select(seqnames, start, end, SNP, A1, A2, BETA, SE, OR, Z, P, N)
      
    } else{
      
      # all other sumstats must have BETA/SE and not OR or Z:
      overlap_df_list[[i]] <-
        gr %>%
        as_tibble() %>%
        dplyr::mutate(
          OR = NA,
          Z = NA,
        ) %>%
        dplyr::select(seqnames, start, end, SNP, A1, A2, BETA, SE, OR, Z, P, N)
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
  
  if (k == 1){
    
    write.table(
      loci,
      here(project_dir,"test_loci", "gwas_filtered.loci"),
      sep = "\t",
      row.names = F,
      quote = F
    )
    
    saveRDS(
      overlap_df,
      file = here(project_dir,"test_loci", "gwas_filtered_loci.rds")
    )
    
  } else if (k == 2){
    
    write.table(
      loci,
      here(project_dir,"test_loci", "gwas_filtered.loci_no_proxies"),
      sep = "\t",
      row.names = F,
      quote = F
    )
    
    saveRDS(
      overlap_df,
      file = here(project_dir,"test_loci", "gwas_filtered_loci_no_proxies.rds")
    )
    
  }
  
}


