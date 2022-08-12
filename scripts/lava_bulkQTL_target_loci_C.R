# LAVA QTL FLD
# Description: run univariate and bivariate tests for eQTL/GWAS traits
# runs the univariate and bivariate tests together for each GWAS/eGene separately, using as target the eqtl in question (arrays).
# tests genome-wide significant genes from QTL sumstats

# Load packages -----------------------------------------------------------

library(here)
library(gtools)
library(LAVA)
library(stringr)
library(data.table)
library(tidyverse)

# Set arguments -----------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# parallelisation chunk
chunk <- c(201:337)
chunk_name <- "201:337"

# list lava sumstats files for eQTLs
eqtls <- list.files(here(project_dir,"GWAS_summary_statistics","LAVA_sumstats","eQTLGen"), recursive = T, full.names = T, pattern = "*.lava.gz")

# remove path and suffix from list of eqtls
eqtl_names <- eqtls %>%
  stringr::str_remove(., here(project_dir, "GWAS_summary_statistics","LAVA_sumstats","eQTLGen","")) %>%
  stringr::str_remove(., ".lava.gz")

# window size flanking the top SNP in the gene
window_size = "100000"

args <-
  list(
    ref_prefix = here(project_dir,"reference_data","g1000_eur","g1000_eur"),
    loc_file = here(project_dir,"test_loci","window_100000", "gene_filtered_bulkQTL.loci"),
    info_file = here(project_dir,"info_files", "input.info.bulkQTL_gwas.txt"),
    sample_overlap_file = here(project_dir, "sample_overlap","bulkQTL_gwas_sample_overlap.txt"),
    gwas_info = here(project_dir, "info_files", "input.info.txt"),
    phenotypes = c(""),
    univar_threshold = c()      # Will be set to 0.05/n_loci
  )

cat("Arguments have been loaded.")

# Load common files -------------------------------------------------------

loci_all <-
  read.loci(args$loc_file) %>%
  dplyr::arrange(LOC)

# "Parallelisation for loci":
loci <-
  loci_all[chunk,]

cat("The number of loci to process is: ", nrow(loci), ". \n")

# create an array of gwas phenotypes to consider in analysis
gwas <-
  read.table(
    args$gwas_info,
    sep = "\t", header = T
  ) %>%
  select("phenotype") %>%
  .[,1]

# Main --------------------------------------------------------------------

# Update univariate threshold
args$univar_threshold <- 0.05/nrow(loci_all)

cat("The univar threshold is: ",args$univar_threshold,". \n")

# create empty list for outputs
univar_one_loop = bivar_one_loop = list()
univar = bivar = list()

# Main loop across all loci (i.e. eGenes)
for (i in 1:nrow(loci)) {
  gene <- loci$LOC[i]
  
  eqtl_gene <- str_c("eQTLGen_",gene)

  if (!eqtl_gene %in% eqtl_names){
    cat("Chunk", i," -- no gene-eqtl pair for", eqtl_gene, "\n")
    next
  }

  cat("Chunk", i, "-- processing locus for", gene,".\n")
  
  # Loop to perform correlations one gwas at a time (to avoid losing SNPs)
  for (j in 1:length(gwas)){
    
    gwas_to_test <- gwas[j]
    
    # Update args (consider one GWAS at a time and the eqtl in question)
    args$phenotypes <-
      c(
        gwas_to_test,
        eqtl_gene
      )
    
    # Load input
    input <-
      try(LAVA::process.input(
        input.info.file = args$info_file,
        sample.overlap.file = args$sample_overlap_file,
        ref.prefix = args$ref_prefix,
        phenos = args$phenotypes
      )
      )
    
    if (inherits(input, "try-error")) {
      next
    }
    
    # Process locus
    locus <-
      LAVA::process.locus(
        loci[i,],
        input
      )
    
    if (!is.null(locus)){
      cat("The locus has been successfully processed for the GWAS: ",gwas_to_test,".\n")
      ls(locus)
    } else{
      
      cat("Chunk ", i, "The locus is null for the GWAS: ",gwas_to_test,".\n")
      next
    }
    
    # It is possible that the locus cannot be defined for various reasons (e.g. too few SNPs),
    # The !is.null(locus) check is necessary before calling the analysis functions.
    if (!is.null(locus)) {
      
      cat("Chunk", i, "-- running uni-bivariate analysis on GWAS ",gwas_to_test," for ", eqtl_gene,".\n")
      
      # extract some general locus info for the output
      loc_info <-
        data.frame(
          locus = locus$id,
          chr = locus$chr,
          start = locus$start,
          stop = locus$stop,
          n_snps = locus$n.snps,
          n_pcs = locus$K
        )
      
      ##### run univ and bivar together and use the eqtl as target
      loc_out <- LAVA::run.univ.bivar(
        locus,
        target = eqtl_gene,
        univ.thresh = args$univar_threshold
      )
      
      # Bind locus info with locus output
      univar_one_loop[[j]] <- loc_info %>%
        bind_cols(loc_out$univ)
      
      if(!is.null(loc_out$bivar)){
      
      bivar_one_loop[[j]] <- loc_info %>%
        bind_cols(loc_out$bivar)
      
      }else{
        
        cat("Chunk ", i, " -- no uni-bivariate test run for GWAS ",gwas_to_test, " on ", eqtl_gene, ".")
        
      }
      
    }
    
  }
      
      # Save locus file into univar and bivar
      univar[[i]] <- bind_rows(univar_one_loop)
      bivar[[i]] <- bind_rows(bivar_one_loop)
      
}

# Save results in single RDS file --------------------------------------------------------------------

saveRDS(
  univar,
  file = here(project_dir,"lava_results","lava_QTLs","univ", str_c("GWAS_with_bulkQTL_all_loci_",chunk_name,".univ.lava.rds"))
)
saveRDS(
  bivar,
  file = here(project_dir,"lava_results","lava_QTLs","bivar", str_c("GWAS_with_bulkQTL_all_loci_",chunk_name,".bivar.lava.rds"))
)

print(str_c(Sys.time(), "Analysis done!"))

