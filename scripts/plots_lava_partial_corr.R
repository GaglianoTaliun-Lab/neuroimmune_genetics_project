library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(qdapTools)
library(purrr)
library(forcats)
library(data.tables)

# Arguments -------------------------------------------------------------------------------------------
project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

source(here(project_dir, "plots_RHR.R"))

pcorr_files <- list.files(here::here(project_dir,"lava_results","lava_GWAS","partial_corr"), pattern = "*partcorr.lava.rds", all.files = T, full.names = F)

locus <- pcorr_files %>% str_extract(., "\\-*\\d+\\.*\\d*")
unique_locus <- unique(locus)

bivar_file <- "AD_schwartzentruber2021:ALS_vanrheenen2016:CD_delange2017:LBD_chia2021:MS_imsgc2019:MSA_scholz2021:PD_nalls2019:SCZ_pardinas2018:UC_delange2017.bivar.lava.rds"

for (nloci in 1:length(unique_locus)) {

pcorr_files_per_locus <- list.files(here::here(project_dir,"lava_results","lava_GWAS","partial_corr"), 
pattern = stringr::str_c("locus_",unique_locus[nloci],"*"), all.files = T, full.names = T)


# Read files --------------------------------------------------------------------------------------------
    bivar_results <-
      readRDS(
        here::here(project_dir,"lava_results","lava_GWAS",bivar_file)
      ) %>%
      purrr::discard(is.null) %>%
      qdapTools::list_df2df() %>%
      dplyr::select(-X1)

# Main loop for plots ------------------------------------------------------------------------------------

    for (i in 1:length(pcorr_files_per_locus)) {
        
        pcorr <- lapply(pcorr_files_per_locus, fread) %>% rbindlist()

        data_to_plot <-
          pcorr %>% 
          dplyr::filter(!is.na(pcor)) %>% 
          dplyr::mutate(locus = as.numeric(unique_locus[nloci])) %>%
          dplyr::select(
            locus, phen1, phen2, pcor, p_cond = p
          ) %>%
          dplyr::inner_join(
            bivar_results %>% 
              dplyr::select(
                locus, contains("phen"), rho, p_uncond = p
              ), 
            by = c("locus", "phen1", "phen2")
          ) %>% 
          tidyr::pivot_longer(
            cols = c("pcor", "p_cond", "rho", "p_uncond")
          ) %>% 
          dplyr::mutate(
            model = 
              case_when(
                name %in% c("pcor", "p_cond") ~ "Conditioned",
                TRUE ~ "Unconditioned"
              ) %>%
              fct_rev(),
            coef = 
              case_when(
                name %in% c("pcor", "rho") ~ "corr",
                TRUE ~ "p"
              )
          ) %>% 
          dplyr::select(-name) %>% 
          tidyr::pivot_wider(
            names_from = c("coef"),
            values_from = c("value")
          )

        plot <- data_to_plot %>% 
          dplyr::mutate(
            rg_fill =
              dplyr::case_when(
                # Different p-value cut-offs depending on conditioning
                p < 0.05 & model == "Conditioned" ~ round(corr, 2),
                p < 0.05/nrow(bivar_results) & model == "Unconditioned" ~ round(corr, 2)
              ),
            p1 = phen1 %>%
              stringr::str_remove(., "_.*"),
            p2 = phen2 %>%
              stringr::str_remove(., "_.*")
          ) %>%
          ggplot2::ggplot(
            ggplot2::aes(
              x = p1,
              y = p2 %>% fct_rev(),
              fill = rg_fill,
              label = round(corr, 2)
            )
          ) +
          ggplot2::geom_tile(colour = "black") +
          ggplot2::geom_text(
            size = 3
          ) +
          ggplot2::facet_grid(rows = vars(model), cols = vars(locus)) +
          ggplot2::coord_fixed() +
          ggplot2::labs(x = "", y = "", fill = "Local genetic correlation (rg)") +
          ggplot2::scale_fill_distiller(
            palette = "RdYlBu",
            limits = c(-1, 1)
          ) +
          ggplot2::scale_colour_manual(values = c("black", "white")) +
          ggplot2::guides(colour = "none") +
          theme_rhr +
          ggplot2::theme(
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank()
          )
    }
       ggsave(here(project_dir,"lava_results","figures",stringr::str_c("partial_corr_locus_",unique_locus[nloci],".pdf")), plot=plot)
}

