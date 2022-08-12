# Description: get gene loci from LD blocks of interest

# Load packages -----------------------------------------------------------

library(here)
library(GenomicRanges)
library(stringr)
library(tidyverse)

# Load data ---------------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# read table provided by LAVA for LD blocks:
ld_blocks <-
  read.table(here(project_dir,"reference_data","blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"),
             sep = " ", header = T)

# read LD blocks used in GWAS LAVA run (v3)
gwas_df <-
  readRDS(here(project_dir,"test_loci","gwas_filtered_loci.v3.rds"))

# read univariate significant results from GWAS (to filter LD blocks)
gwas_univ_sign <- read.table(
  here(project_dir,"lava_results","tsv_univ_bivar_GWAS","phenotypes_significant_results_univ.v3.tsv"),
  sep = "\t", header = T)

# can obtain two window sizes
window_size <- c(50000, 100000)

# load gene table reference GRCh37
ref <- rtracklayer::import(here(project_dir,"reference_data","Homo_sapiens.GRCh37.87.gtf"))
ref <- ref %>% GenomeInfoDb::keepSeqlevels(c(1:22), pruning.mode = "coarse")
ref <- ref[ref$type == "gene"]

# Main --------------------------------------------------------------------

# Add locus ID and rename remaining columns to fit LAVA requirements; retain only LD blocks in GWAS LAVA run (v3)
ld_blocks <-
  ld_blocks %>%
  dplyr::rename_with(
    .fn = stringr::str_to_upper,
    .cols = everything()
  ) %>%
  dplyr::mutate(
    LOC = dplyr::row_number()
  ) %>%
  dplyr::select(LOC, everything()) %>%
  dplyr::filter(LOC %in% unique(gwas_univ_sign$locus))

# Convert ld_blocks to granges
ld_blocks_gr <-
  ld_blocks %>%
  GenomicRanges::makeGRangesFromDataFrame(
    .,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqinfo = NULL,
    seqnames.field = "CHR",
    start.field = "START",
    end.field = "STOP"
  )

# Overlap granges objects: the gene reference table and LD blocks filtered
overlap <-
  GenomicRanges::findOverlaps(ld_blocks_gr, ref) %>%
  tibble::as_tibble()

# Generate df of genes that overlap which LD blocks
overlap_df <-
  window_size %>%
  lapply(., function(window){
    
    tibble::tibble(
      locus = ld_blocks_gr[overlap$queryHits]$LOC,
      gene_id = ref[overlap$subjectHits]$gene_id,
      gene_name =
        ref[overlap$subjectHits]$gene_name %>%
        stringr::str_replace_all("-", "_") %>%
        stringr::str_replace_all("\\.", "_"),
      gene_biotype = ref[overlap$subjectHits]$gene_biotype %>% as.factor(),
      chr = ref[overlap$subjectHits] %>% GenomeInfoDb::seqnames() %>% as.character(),
      strand = ref[overlap$subjectHits] %>% BiocGenerics::strand() %>% as.character(),
      start = c(ref[overlap$subjectHits] %>% BiocGenerics::start()),
      end = ref[overlap$subjectHits] %>% BiocGenerics::end(),
      width = ref[overlap$subjectHits] %>% BiocGenerics::width()
    ) %>%
      # If there are two versions of a gene, choose the longest
      dplyr::group_by(gene_name) %>%
      dplyr::top_n(1, wt = width) %>%
      dplyr::filter(
        gene_biotype %in% c("protein_coding", "antisense", "lincRNA")
      ) %>%
      # Add window
      dplyr::mutate(
        start_window =
          case_when(
            start - window <= 0 ~ 0,
            TRUE ~ start - window
          ),
        end_window = end + window
      )
    
  })

names(overlap_df) <- str_c("window_", as.integer(window_size))

# Create locus file 
loci <-
  overlap_df %>%
  lapply(., function(df){
    
    df %>%
      dplyr::ungroup() %>%
      dplyr::distinct(
        gene_id, chr, start_window, end_window
      ) %>%
      dplyr::select(
        LOC = gene_id,
        CHR = chr,
        START = start_window,
        STOP = end_window
      )
    
  })


# Generate df of gwas loci that overlap genes/ld blocks from QTLs

# create granges object from GWAS loci
gwas_gr <-
  gwas_df %>%
  GenomicRanges::makeGRangesFromDataFrame(
    .,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqinfo = NULL,
    seqnames.field = "GWAS_CHR",
    start.field = "GWAS_start",
    end.field = "GWAS_end"
  )

# create granges object from gene LD blocks and reference gene table
locus_gr <-
  overlap_df %>%
  lapply(., function(df){
    
    df %>%
      dplyr::ungroup() %>%
      dplyr::distinct(
        locus, gene_id, chr, start_window, end_window
      ) %>%
      GenomicRanges::makeGRangesFromDataFrame(
        .,
        keep.extra.columns = TRUE,
        ignore.strand = TRUE,
        seqinfo = NULL,
        seqnames.field = "chr",
        start.field = "start_window",
        end.field = "end_window"
      )
    
  })

# get the overlap of both granges objects (locus_gr and gwas_gr)
overlap_gwas_gene <-
  locus_gr %>%
  lapply(., function(gr){
    
    gr %>%
      GenomicRanges::findOverlaps(gwas_gr, .) %>%
      tibble::as_tibble()
    
  })

# create empty list for output results (data frame of gene-gwas overlap)
gwas_gene_overlap_df <-
  setNames(
    vector(length = length(overlap_gwas_gene), mode = "list"),
    nm = names(overlap_gwas_gene)
  )

# fill in list with the gene-gwas overlaps
for(window in names(overlap_gwas_gene)){
  
  gwas_gene_overlap_df[[window]] <-
    gwas_df %>%
    dplyr::slice(overlap_gwas_gene[[window]]$queryHits) %>%
    dplyr::select(contains("GWAS")) %>%
    dplyr::bind_cols(
      locus_gr[[window]] %>%
        as_tibble() %>%
        dplyr::rename_with(
          ~ stringr::str_c("GENE", .x, sep = "_")
        ) %>%
        dplyr::slice(overlap_gwas_gene[[window]]$subjectHits)
    ) %>%
    dplyr::select(contains("GENE"), contains("GWAS"), -contains("strand")) %>%
    dplyr::rename_with(
      ~ stringr::str_replace(.x,
                             pattern = "seqnames",
                             replacement = "CHR"),
      .col = dplyr::ends_with("seqnames")
    )
  
}


# Save data ---------------------------------------------------------------

out_dir <- file.path("/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project/test_loci", names(loci))

for(path in out_dir){
  
  file <- basename(path)
  
  write_delim(
    loci[[file]],
    delim = "\t",
    file = file.path(path, "gene_filtered.loci")
  )
  saveRDS(
    overlap_df[[file]],
    file = file.path(path, "gene_filtered_loci.rds")
  )
  
  saveRDS(
    gwas_gene_overlap_df[[file]],
    file = file.path(path, "gene_gwas_filtered_loci.rds")
  )
  
}

