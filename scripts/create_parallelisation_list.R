# Description: short script to create list of loci for LAVA parallelisation (arrays)

library(tidyverse)
library(stringr)
library(here)

project_dir="/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

start=1
end=2553 # total number of unique genes across all cell types
nplus=100

paral_list <- data.frame(chunk_start = seq(start,end, nplus), chunk_end = seq(100, end+100, nplus))
paral_list[nrow(paral_list),2] = end
paral_list <- paral_list %>%
  mutate(., chunk_name = str_c(chunk_start, ":", chunk_end)) %>%
  select(
    chunk_name,
    chunk_start,
    chunk_end
  )

write.table(paral_list, here(project_dir,"parallelisation_list.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
