library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(data.table)

# Arguments ----------------------------------------------------------------------------------------------------------------------
project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

source(here(project_dir, "plots_RHR.R"))

multireg_files <- list.files(here(project_dir,"lava_results","lava_GWAS","multi_reg"), pattern = "*multireg.lava.rds", all.files = T, full.names = F)

locus <- multireg_files %>% str_extract(., ".*target") %>% str_remove(., "_target") %>% str_replace(., "_", " ") %>% str_remove(., "locus ")
unique_locus <- unique(locus)

# Main ------------------------------------------------------------------------------------------------------------------------

# Loop across unique LAVA loci
for (nloci in 1:length(unique_locus)) {
  
  # list multiple regression results paths
  multireg_files_full <- list.files(here::here(project_dir, "lava_results", "lava_GWAS", "multi_reg"),
                                         pattern = stringr::str_c("locus_", unique_locus[nloci], "*"), all.files = T, full.names = T)
  

  # create empty list with length equal to the number of multiple regression results files
  data_to_plot <- vector(mode = "list", length = length(multireg_files_full))
  
  # Loop across each multiple regression file
  for (i in 1:length(multireg_files_full)) {
    
    # list multiple regression files
    multireg_file <- list.files(here::here(project_dir, "lava_results", "lava_GWAS", "multi_reg"),
                               pattern = stringr::str_c("locus_", unique_locus[nloci], "*"), all.files = T, full.names = F) %>% .[i]
    
    # read multiple regression results in list
    data_to_plot[[i]] <- fread(multireg_files_full[i])
    

    # Define outcome, predictors and models for each multiple regression to use in plots
    outcome <- multireg_file %>% str_remove_all(., ".*target_") %>% str_extract(., ".*predictors") %>% str_remove(., "_predictors")
    predictors_out <- multireg_file %>% str_remove_all(., ".*target_") %>% str_extract(., "predictors_.*") %>% str_remove(., ".multireg.lava.rds") %>% str_remove(., "predictors_") 
    predictors <- predictors_out %>% str_replace(., "_", " + ")
    model <- multireg_file %>% str_remove_all(., ".*target_") %>% str_remove(., ".multireg.lava.rds") %>% str_replace(., "_predictors_", " ~ ") %>% str_replace(., "_", " + ")
    
    # include new column with significance symbols for plots
    data_to_plot[[i]] <-   
      data_to_plot[[i]] %>%  
      dplyr::mutate(
        p_text = case_when(
          p < 0.001 ~ "***",
          p < 0.01 ~ "**",
          p < 0.05 ~ "*",
          p >= 0.05 ~ "ns"
        ),
        locus = as.numeric(unique_locus[nloci]),
        model = str_c(model,""),
        predictors = stringr::str_remove(predictors, "_.*"), outcome = stringr::str_remove(outcome, "_.*")
      )
    
    # wrangle data for plotting
    data_to_plot[[i]] <- 
      data_to_plot[[i]] %>% 
      dplyr::inner_join(
        data_to_plot[[i]] %>% 
          dplyr::group_by(outcome, locus) %>% 
          dplyr::select(outcome, locus, model)
      )
    
  }
  
  all_data_plot <- data.table::rbindlist(data_to_plot)
  
  # Main plots code --------------------------------------------------------------------------------------------------
    a <- 
      all_data_plot %>% 
      dplyr::mutate(
        locus_model =
          str_c(locus, model, sep = "\n")
      ) %>%  
      dplyr::distinct(
        locus, locus_model, r2, r2.lower, r2.upper
      ) %>% 
      ggplot2::ggplot(
        ggplot2::aes(
          x = locus_model,
          y = r2,
          ymin = r2.lower,
          ymax = r2.upper
        )
      ) +
      geom_col() + 
      geom_errorbar(width=.2) +
      labs(
        x = "Locus and model",
        y = "Multivariate r2"
      ) +
      theme_rhr
    
    b <- 
      all_data_plot %>% 
      ggplot2::ggplot(
        ggplot2::aes(
          x = predictors, 
          y = gamma, 
          ymin = gamma.lower, 
          ymax = gamma.upper
        )
      ) +
      ggplot2::geom_pointrange(
        colour = "#888888", fatten = 0.25
      ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::geom_text(
        aes(
          y = gamma.upper + (0.1 * gamma.upper), 
          label = p_text
        )
      ) +
      ggplot2::facet_wrap(
        vars(locus, model), 
        scales = "free_x",
        nrow = 1
      ) +
      ggplot2::labs(x = "Predictors", y = "Standardised multiple\nregression coefficient") +
      theme_rhr
    

    # plot a and b ggplot elements into one figure
    plot <- 
      cowplot::plot_grid(
        b,a, 
        ncol = 1, 
        axis = "lr", 
        align = "v"
      )
  
  # save plots  
  ggsave(here(project_dir,"lava_results","figures",stringr::str_c("multi_reg_locus_",unique_locus[nloci],".pdf")), plot=plot)

}
