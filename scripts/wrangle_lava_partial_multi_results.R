library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

pcorr <-
  list.files(
    here::here(
      "results", 
      "03_partial_corr_multi_reg"
      ),
    pattern = ".partcorr.lava.rds", 
    full.names = T
  ) %>% 
  readRDS() %>% 
  lapply(., function(list){
    
    list %>% 
      qdapTools::list_df2df(col1 = "list_name")
    
  }) %>% 
  qdapTools::list_df2df(col1 = "list_name_2") %>% 
  dplyr::select(-contains("list_name")) %>% 
  as_tibble()
  
multireg <-
  list.files(
    here::here(
      "results", 
      "03_partial_corr_multi_reg"
      ),
    pattern = ".multireg.lava.rds", 
    full.names = T
  ) %>% 
  readRDS() %>% 
  lapply(., function(x){
    
    x %>% 
      lapply(., function(y){
        
        y %>% 
          lapply(., function(z){
            
            z %>% 
              qdapTools::list_df2df(col1 = "model_number")
            
          }) %>% 
          qdapTools::list_df2df(col1 = "list_name")
        
      }) %>% 
      qdapTools::list_df2df()
    
  }) %>% 
    qdapTools::list_df2df(col1 = "locus") %>% 
    dplyr::select(-X1) %>% 
    dplyr::mutate(
      model_type = 
        case_when(
          locus == 1719 & list_name == "L1" ~ "intermediate",
          locus == 1719 & list_name == "L2" ~ "full",
          TRUE ~ "full"
          )
    ) %>% 
    dplyr::select(
      locus, contains("model"), 
      outcome, predictors, 
      everything(), -contains("list")
      )
