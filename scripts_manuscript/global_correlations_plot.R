# Description: script to obtain heat map plots of the global genetic correlations across traits.

# Libraries ----------------------------------------------------------------------------------------

library(here)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gghighlight)

# Arguments ----------------------------------------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

correlations <- read.table(here(project_dir,"ldsc_corr","ldsc_correlations.txt", sep = "\t", header = T))

# Wrangle data ----------------------------------------------------------------------------------------

# round rg to 2 digits
correlations <- correlations %>%
  mutate(rg = round(rg, digits = 2))

# remove unwanted phenotypes, with proxies:
correlations_prox <- correlations %>% 
  filter(p1 != "AD_wightman2021") %>%
  filter(p2 != "AD_wightman2021") %>%
  filter(p1 != "AD_kunkle2019") %>%
  filter(p2 != "AD_kunkle2019") %>%
  filter(p1 != "SCZ_lam2019") %>%
  filter(p2 != "SCZ_lam2019") %>%
  filter(p1 != "MS_sawcer2011") %>%
  filter(p2 != "MS_sawcer2011") %>%
  filter(p1 != "PD_blauwendraat2019") %>%
  filter(p2 != "PD_blauwendraat2019") %>%
  filter(p1 != "ALS_vanrheenen2016") %>%
  filter(p2 != "ALS_vanrheenen2016") %>%
  mutate(p1 = str_remove(p1, "_.*"), p2 = str_remove(p2, "_.*"))

# remove unwanted phenotypes, without proxies:
correlations_noprox <- correlations %>% 
  filter(p1 != "AD_wightman2021") %>%
  filter(p2 != "AD_wightman2021") %>%
  filter(p1 != "AD_schwartzentruber2021") %>%
  filter(p2 != "AD_schwartzentruber2021") %>%
  filter(p1 != "SCZ_lam2019") %>%
  filter(p2 != "SCZ_lam2019") %>%
  filter(p1 != "MS_sawcer2011") %>%
  filter(p2 != "MS_sawcer2011") %>%
  filter(p1 != "PD_nalls2019") %>%
  filter(p2 != "PD_nalls2019") %>%
  filter(p1 != "ALS_vanrheenen2016") %>%
  filter(p2 != "ALS_vanrheenen2016") %>%
  filter(p1 != "NA") %>%
  mutate(p1 = str_remove(p1, "_.*"), p2 = str_remove(p2, "_.*"))

correlations_list <- list(correlations_prox, correlations_noprox)

# Main loop across traits ----------------------------------------------------------------------------------------

for (i in 1:length(correlations_list)) {
  
  # set p-value threshold
  n = unique(correlations_list[[i]]$p1) %>% length(.)
  n_tests = (n*(n-1)/2)
  p_threshold <- 0.05/n_tests
  
  # save table with significant correlations
  corr_sign <- filter(correlations_list[[i]], p <= p_threshold)
  # write.table(corr_sign, paste0(out_dir,"significant_global_correlations.tsv"), sep = "\t", quote = F, row.names = F)
  
  # factor to order phenotypes:
  correlations_list[[i]]$p1.ord <- factor(correlations_list[[i]]$p1, levels= rev(c("CD","UC","MS","SCZ","ALS","FTD","MSA","PD","LBD","AD")))
  correlations_list[[i]]$p2.ord <- factor(correlations_list[[i]]$p2, levels= c("CD","UC","MS","SCZ","ALS","FTD","MSA","PD","LBD","AD"))
  
  # create matrix from long data
  correlations_matrix <- correlations_list[[i]][,c(13,14,3)] %>% spread(., p2.ord, rg)
  
  # create plot with matrix and corrplot library
  # currently not working because it shows funny half matrix
  # corrplot(correlations_matrix, type = "lower", 
  #         tl.col = "black", tl.srt = 45)
  # http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
  
  # create plot
  ggplot(data = correlations_list[[i]], mapping = aes(x = p1.ord, y = p2.ord, fill = rg)) +
    geom_tile(colour = "black", lwd = 0.5, linetype = 1) +
    xlab(label = "") +
    ylab(label = "") +
    theme_classic() +
    scale_fill_gradient2(low = "purple",
                         high = "red3", na.value = "white",
                         limits = c(-1, 1.08)) +
    guides(fill = guide_colourbar(title = "Global Correlations (rg)"), x = guide_axis(angle = 90)) + 
    theme(axis.text.x = element_text(size = 16, angle = 90, hjust = .5, vjust = .5),
          axis.text.y = element_text(size = 16, hjust = .5, vjust = .5),
          axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
    coord_fixed() +
    geom_text(aes(label=rg), size = 6) +
    gghighlight(p < p_threshold,
                unhighlighted_params = list(colour = alpha("white")))
  
  if (i == 1) {
    
    ggsave(here(project_dir,"ldsc_corr","figures","global_correlations_significant_with_proxies.pdf"), width = 40, height = 30, units = "cm")
    
  }else if (i == 2) {
    
    ggsave(here(project_dir,"ldsc_corr","figures","global_correlations_significant_no_proxies.pdf"), width = 40, height = 30, units = "cm")
    
  }
  
}
