# Description: script to wrangle results of LAVA with ONEK1K QTLs

# Packages ----------------------------------------------------------------

library(here)
library(purrr)
library(dplyr)
library(tidyr)
library(stringr)
library(qdapTools)
library(qvalue)

# Arguments ----------------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

date=str_remove_all(Sys.Date(), "-")

results <- 
  setNames(
    vector(mode = "list", length = 2),
    nm = c("univ", "bivar")
  ) %>% 
  lapply(., function(x){
    
    setNames(
      vector(mode = "list", length = 1),
      nm = "window_100000"
    )
    
  })

window = "window_100000"

# Main ----------------------------------------------------------------

for(test in names(results)){
  
  gene_filtered_loci <-
    readRDS(
      here(project_dir, "test_loci", window, "onek1k_filtered_loci.rds")
    )
  
  results[[test]][[window]] <- 
    setNames(
      object = 
        list.files(
          path = 
            here(project_dir,"lava_results","lava_onek1k_QTLs",test),
          pattern = str_c("PD_replication.",test,".lava.rds"), 
          full.names = T
        ) %>% 
        lapply(., function(file) readRDS(file)),
      nm = 
        list.files(
          path = 
            here(project_dir,"lava_results","lava_onek1k_QTLs",test),
          pattern = str_c("PD_replication.",test,".lava.rds")
        ) %>% 
        str_remove(., ".lava.rds") %>%
        str_remove(., here(project_dir,"lava_results","lava_onek1k_QTLs",test,"")) %>%
        str_remove(., "GWAS_QTLs_") %>%
        str_remove(., "_PD_replication")
    ) %>% 
    purrr::discard(is.null)
  
  if (test == "univ") {
    
    results_wrangled <- 
      results[[test]][[window]] %>%
      discard(is.null) %>%
      lapply(., function(x) bind_rows(x))
    
    for (i in 1:length(results_wrangled)){
      results_wrangled[[i]] <-
        results_wrangled[[i]] %>%
        mutate(., qtl_dataset = names(results_wrangled[i]) %>% stringr::str_remove(., ".univ"))
    }
    
    results_wrangled <-
      results_wrangled %>%
      bind_rows(.) %>%
      dplyr::rename(
        gene_locus = locus
      ) %>%
      dplyr::inner_join(
        gene_filtered_loci %>% 
          dplyr::select(
            ld_block = locus, gene_id, gene_name
          ),
        by = c("gene_locus" = "gene_id")
      ) %>% 
      dplyr::select(
        gene_locus, gene_name, everything()
      ) %>% 
      dplyr::arrange(chr, start) %>% 
      select(!ld_block)	%>%
      as_tibble() %>% unique()
    
    write.table(results_wrangled, here(project_dir,"lava_results","tsv_univ_bivar_QTLs",stringr::str_c("onek1k_PD_replication_",test,"_",window,"_",date,".tsv")),
                sep = "\t", quote=F, row.names=F)
    
  } else if (test == "bivar") {
    
    results_wrangled <- 
      results[[test]][[window]] %>%
      discard(is.null) %>%
      bind_rows(.) %>%
      dplyr::rename(
        gene_locus = locus
      ) %>%
      dplyr::inner_join(
        gene_filtered_loci %>% 
          dplyr::select(
            ld_block = locus, gene_id, gene_name
          ),
        by = c("gene_locus" = "gene_id")
      ) %>% 
      dplyr::select(
        gene_locus, gene_name, everything()
      ) %>% 
      dplyr::arrange(chr, start) %>%
      select(!ld_block) %>% 
      as_tibble() %>% unique()
    
    write.table(results_wrangled, here(project_dir,"lava_results","tsv_univ_bivar_QTLs",stringr::str_c("onek1k_PD_replication_",test,"_",window,"_",date,".tsv")),
                sep = "\t", quote=F, row.names=F)
    
    # Bonferroni correction
    n_tests <- nrow(results_wrangled)
    pvalue_bonf <- 0.05/n_tests
    results_wrangled %>% 
      dplyr::filter(., p <= pvalue_bonf) %>%
      write.table(., here(project_dir,"lava_results","tsv_univ_bivar_QTLs",stringr::str_c("bonferroni_significant_onek1k_PD_replication_",test,"_",window,"_",date,".tsv")),
                  sep = "\t", quote=F, row.names=F)
    
    print(str_c("Number of total bivariate tests: ", n_tests, ". The bonferroni corrected p-value = ", pvalue_bonf, "."))
    
  }
  
}
