# Description: script to obtain the number of unique genes across all treatments from the sc-eQTLs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(devtools)
library(qdapTools)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(here)

# devtools::install_github("RHReynolds/colochelpR")
# devtools::install_github("RHReynolds/rutils")
# library(colochelpR)
# library(rutils)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Arguments
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"
lava_ext = ".lava.gz"

# empty df:
all_unique <- data.frame()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (type in c("B","CD4T","CD8T","monocyte","NK")){

  qtls_list <-
    setNames(
      object = 
        list.files(here(project_dir,"GWAS_summary_statistics","1M-scBloodNL","treatments_eqtls"),
                   pattern = stringr::str_c(type,"_eQTLs.txt.gz"), full.names = T) %>%
        lapply(., function(qtl)fread(qtl, fill=TRUE)
               ),
      nm =
        list.files(here(project_dir,"GWAS_summary_statistics","1M-scBloodNL","treatments_eqtls"),
                   pattern = stringr::str_c(type,"_eQTLs.txt.gz"), full.names = F) %>%
        stringr::str_remove(., "_eQTLs.txt.gz")
      )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  qtls_list <-
    qtls_list %>%
    lapply(., function(qtl){
        filter(qtl, PValue < 5e-08) %>%
        select(
          PValue,
          ProbeName
        )
    })

  unique_genes <- data.table::rbindlist(qtls_list) %>%
    distinct(., ProbeName, .keep_all = TRUE)

  cat("The number of unique genes for ",type," eQTLs is: ",nrow(unique_genes),".")

  all_unique <- rbind(all_unique, unique_genes)

}

# output number of unique genes across all treatments:
all_unique_final <- all_unique %>%
  distinct(., ProbeName, .keep_all = TRUE)

cat("The number of unique genes across all treatments is: ",nrow(all_unique_final),".")

