# Description: generate scatter plots to look at the correlations between ALS (2016,2021) GWAS sumstats
# and SCZ (2018).

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ALS2021 <- read.table(here(project_dir,"GWAS_summary_statistics","ldsc_sumstats","ALS_vanrheenen2021.sumstats.gz"),
  sep = "\t", header = T, fill = TRUE) %>%
  rename(Z_ALS2021 = Z)

ALS2016	<- read.table(here(project_dir,"GWAS_summary_statistics","ldsc_sumstats","ALS_vanrheenen2016.sumstats.gz"),
  sep = "\t", header = T, fill = TRUE) %>%
  rename(Z_ALS2016 = Z)

SCZ2018	<- read.table(here(project_dir,"GWAS_summary_statistics","ldsc_sumstats","SCZ_pardinas2018.sumstats.gz"),
  sep = "\t", header = T, fill = TRUE) %>%
  rename(Z_SCZ2018 = Z)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ALS2021_ALS2016 <- inner_join(ALS2021, ALS2016, by = "SNP")
ALS2021_SCZ2018 <- inner_join(ALS2021, SCZ2018, by = "SNP")
ALS2016_SCZ2018 <- inner_join(ALS2016, SCZ2018, by = "SNP")

a <- ALS2021_ALS2016 %>% ggplot(aes(x = Z_ALS2021, y = Z_ALS2016)) +
  geom_point()

b <- ALS2021_SCZ2018 %>% ggplot(aes(x = Z_ALS2021, y = Z_SCZ2018)) +
  geom_point()

c <- ALS2016_SCZ2018 %>% ggplot(aes(x = Z_ALS2016, y = Z_SCZ2018)) +
  geom_point()

# plot a,b,c ggplot elements into one figure
plot <- 
  cowplot::plot_grid(
    a,b,c, 
    ncol = 3, 
    axis = "lr", 
    align = "v"
  )
  
# save plot
ggsave(here(project_dir,"ldsc_corr","figures","scatter_plot_GWAS_sumstats_ALS_SCZ.pdf"), plot=plot)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Correlation test
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cor_ALS2016_ALS2021 <- cor.test(ALS2021_ALS2016$Z_ALS2021, ALS2021_ALS2016$Z_ALS2016)
cor_ALS2021_SCZ2018 <- cor.test(ALS2021_SCZ2018$Z_ALS2021, ALS2021_SCZ2018$Z_SCZ2018)
cor_ALS2016_SCZ_2018 <- cor.test(ALS2016_SCZ2018$Z_ALS2016, ALS2016_SCZ2018$Z_SCZ2018)

cat("Correlation between ALS 2016 and ALS 2021 is: ",cor_ALS2016_ALS2021$estimate[[1]]," with p-value = ",cor_ALS2016_ALS2021$p.value,". \n")
cat("Correlation between SCZ 2018 and ALS 2021 is: ",cor_ALS2021_SCZ2018$estimate[[1]]," with p-value = ",cor_ALS2021_SCZ2018$p.value,". \n")
cat("Correlation between SCZ 2018 and ALS 2016 is: ",cor_ALS2016_SCZ_2018$estimate[[1]]," with p-value = ",cor_ALS2016_SCZ_2018$p.value,". \n")
