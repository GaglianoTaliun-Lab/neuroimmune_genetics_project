# Description: get gene loci from LD blocks of interest
# This script is different from 'get_genic_regions.R' in that here it includes an extra filter for
# keeping only genes in the output files that were genome-wide significant in the QTL lava.gz sumstats.

# Load packages -----------------------------------------------------------

library(here)
library(GenomicRanges)
library(stringr)
library(data.table)
library(tidyverse)

# Load data ---------------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# read table provided by LAVA for LD blocks:
ld_blocks <-
  read.table(here(project_dir,"reference_data","blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"),
             sep = " ", header = T)

# read LD blocks used in GWAS LAVA run (v3)
# gwas_df <-
#   readRDS(here(project_dir,"test_loci","gwas_filtered_loci.v3.rds"))

# window size flanking each gene = 100kb
window <- 100000

# load gene table reference GRCh37
ref <- rtracklayer::import(here(project_dir,"reference_data","Homo_sapiens.GRCh37.87.gtf"))
ref <- ref %>% GenomeInfoDb::keepSeqlevels(c(1:22), pruning.mode = "coarse")
ref <- ref[ref$type == "gene"]

# Main --------------------------------------------------------------------

# read QTL sumstats and create a list of unique genes that passed the genome-wide threshold (5e-08)
qtl_list <-
  setNames(
    object =
      list.files(
        path = here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","1MscBloodNL"),
        pattern = "*.lava.gz",
        recursive = T,
        full.names = T
      ) %>%
      lapply(., function(file) fread(file)),
    nm =
      list.files(
        path = here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","QTLs_all"),
        pattern = "*.lava.gz",
        recursive = T,
        full.names = F
      ) %>%
      str_remove(., ".lava.gz")
  )

# filter for genes with genome-wide significant p-values
qtl_sign <- lapply(qtl_list, function(sumstats) filter(sumstats, P <= 5e-08))

unique_genes <- bind_rows(qtl_sign) %>%
  distinct(., GENE)

cat("The number of unique genes across all QTL genome-wide significant signals is: ", nrow(unique_genes), ". \n")

# save table of unique genes:
write.table(unique_genes, here(project_dir, "test_loci", "window_100000","unique_GWS_genes.txt"), sep = "\t", row.names = F, quote = F)

# Add locus ID and rename remaining columns to fit LAVA requirements
  ld_blocks %>%
  dplyr::rename_with(
    .fn = stringr::str_to_upper,
    .cols = everything()
  ) %>%
  dplyr::mutate(
    LOC = dplyr::row_number()
  ) %>%
  dplyr::select(LOC, everything())

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
  # Filter for genome-wide significant QTL genes
  dplyr::filter(
    gene_id %in% unique_genes$GENE
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

# Create locus file 
loci <-
  overlap_df %>%
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

# Generate df of gwas loci that overlap genes/ld blocks from QTLs

# create granges object from GWAS loci
# gwas_gr <-
#   gwas_df %>%
#   GenomicRanges::makeGRangesFromDataFrame(
#     .,
#     keep.extra.columns = TRUE,
#     ignore.strand = TRUE,
#     seqinfo = NULL,
#     seqnames.field = "GWAS_CHR",
#     start.field = "GWAS_start",
#     end.field = "GWAS_end"
#   )

# create granges object from gene LD blocks and reference gene table
# locus_gr <-
#   overlap_df %>%
#   dplyr::ungroup() %>%
#   dplyr::distinct(
#     locus, gene_id, chr, start_window, end_window
#   ) %>%
#   GenomicRanges::makeGRangesFromDataFrame(
#     .,
#     keep.extra.columns = TRUE,
#     ignore.strand = TRUE,
#     seqinfo = NULL,
#     seqnames.field = "chr",
#     start.field = "start_window",
#     end.field = "end_window"
#   )

# get the overlap of both granges objects (locus_gr and gwas_gr)
# overlap_gwas_gene <-
#   locus_gr %>%
#   GenomicRanges::findOverlaps(gwas_gr, .) %>%
#   tibble::as_tibble()

# fill in list with the gene-gwas overlaps
# gwas_gene_overlap_df  <-
#   gwas_df %>%
#   dplyr::slice(overlap_gwas_gene$queryHits) %>%
#   dplyr::select(contains("GWAS")) %>%
#   dplyr::bind_cols(
#     locus_gr %>%
#       as_tibble() %>%
#       dplyr::rename_with(
#         ~ stringr::str_c("GENE", .x, sep = "_")
#       ) %>%
#       dplyr::slice(overlap_gwas_gene$subjectHits)
#   ) %>%
#   dplyr::select(contains("GENE"), contains("GWAS"), -contains("strand")) %>%
#   dplyr::rename_with(
#     ~ stringr::str_replace(.x,
#                            pattern = "seqnames",
#                            replacement = "CHR"),
#     .col = dplyr::ends_with("seqnames")
#   ) %>%
#   dplyr::filter(
#     GENE_gene_id %in% unique_genes$GENE
#   )

# Save data ---------------------------------------------------------------

write_delim(
  loci,
  delim = "\t",
  file = here(project_dir, "test_loci", "window_100000", "gene_filtered.loci")
)

saveRDS(
  overlap_df,
  file = here(project_dir, "test_loci", "window_100000", "gene_filtered_loci.rds")
)

# saveRDS(
#   gwas_gene_overlap_df,
#   file = here(project_dir, "test_loci", "window_100000", "gene_gwas_filtered_loci.rds")
# )


