#' Plot global correlations
#'
#' @description This function will plot global genetic correlations between
#'   phenotypes (as determined using LDSC) as a lower triangle heatmap.
#'   Significant correlations will be determined with multiple test corrections
#'   applied (Bonferroni). Number of tests = to number of unique combinations
#'   between phenotypes (excluding comparisons between the same phenotype).
#'
#' @param global_corr a `data.frame` or [tibble][tibble::tbl_df-class] object,
#'   with the following columns (all of which are output by LDSC genetic
#'   correlation analyses):
#'  \itemize{
#'  \item `p1`: name of phenotype 1
#'  \item `p2`: name of phenotype 2
#'  \item `rg`: the estimated genetic correlation
#'  \item `p`: p-value of the genetic correlation
#'  }
#' @param n_phenotypes `integer` vector indicating number of phenotypes run
#'
#' @return `ggplot` displaying the genetic correlations between phenotypes.
#' \itemize{
#' \item x and y-axis display the phenotypes.
#' \item Significant negative and positive correlations are indicated by blue
#' and red fill, respectively. Non-significant correlations (p >= 0.05/n_tests)
#' have a grey fill.
#' }
#'
#' @export
#'
#' @references \itemize{ \item Bulik-Sullivan et al. (2015) An atlas of genetic
#'   correlations across human diseases and traits \emph{Nature Genetics}, 2015
#'   Nov;47(11):1236-41. \url{https://www.nature.com/articles/ng.3406} PMID:
#'   26414676 }

plot_global_corr <-
  function(
    global_corr,
    n_phenotypes
  ) {
    
    # Only need first half of matrix, thus must extract appropriate rows from dataframe
    n <- n_phenotypes
    
    for (i in 1:n) {
      
      # General formula for extracting appropriate indices
      index <- (i * n - (n - i)):(i * n)
      
      if (i == 1) {
        indices <- index
      } else {
        indices <- c(indices, index)
      }
    }
    
    # Determine number of combinations for multiple test correction
    n_combn <-
      length(indices) - n_phenotypes
    
    global_corr %>%
      dplyr::mutate(
        rg_fill =
          dplyr::case_when(
            p < 0.05 / n_combn ~ round(rg, 2)
          ),
        p1 = p1 %>%
          stringr::str_replace_all("[:digit:]", "") %>%
          stringr::str_remove("\\..*"),
        p2 = p2 %>%
          stringr::str_replace_all("[:digit:]", "") %>%
          stringr::str_remove("\\..*")
      ) %>%
      dplyr::slice(indices) %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = ord_p1,
          y = ord_p2 %>% fct_rev(),
          fill = rg_fill,
          label = round(rg, 2)
        )
      ) +
      ggplot2::geom_tile(colour = "black") +
      ggplot2::geom_text(
        size = 8
      ) +
      ggplot2::coord_fixed() +
      ggplot2::labs(x = "", y = "", fill = "Correlation (rg)") +
      ggplot2::scale_fill_gradient2(low = "purple",
                                    high = "red3", na.value = "grey",
                                    limits = c(-1, 1.08)) +
      ggplot2::scale_colour_manual(values = c("black", "white")) +
      ggplot2::guides(colour = "none") +
      theme_classic() +
      ggplot2::theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        # axis.ticks = element_blank(),
        axis.text.x = element_text(size = 22, hjust = .5, vjust = .5),
        axis.text.y = element_text(size = 22, hjust = .5, vjust = .5),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 22)
      )
  }
