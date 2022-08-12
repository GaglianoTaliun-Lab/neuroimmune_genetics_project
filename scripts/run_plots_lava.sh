#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=0:30:00
#SBATCH --job-name=plots-lava
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=3G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

# arguments:
# arg 1-x = phenotypes that were included in the LAVA analysis

Rscript plots_lava.R AD_schwartzentruber2021 ALS_vanrheenen2021 CD_delange2017 FTD_ferrari2014 LBD_chia2021 MS_imsgc2019 PD_nalls2019 SCZ_pardinas2018 UC_delange2017
# Rscript plots_lava_no_proxies.R AD_kunkle2019 ALS_vanrheenen2021 CD_delange2017 FTD_ferrari2014 LBD_chia2021 MS_imsgc2019 PD_blauwendraat2019 SCZ_pardinas2018 UC_delange2017

# No arguments for the following scripts:
# Rscript plots_lava_multi_reg.R
# Rscript plots_lava_partial_corr.R
# Rscript plots_lava_qtls.R
