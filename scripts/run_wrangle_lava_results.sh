#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=0:30:00
#SBATCH --job-name=wrangle-lava
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

# arguments:
# arg 1 = prefix of bivar or univ test lava output without path

# Rscript wrangle_lava_results.R AD_schwartzentruber2021:ALS_vanrheenen2021:CD_delange2017:FTD_ferrari2014:LBD_chia2021:MS_imsgc2019:PD_nalls2019:SCZ_pardinas2018:UC_delange2017
Rscript wrangle_lava_results_no_proxies.R AD_kunkle2019:ALS_vanrheenen2021:CD_delange2017:FTD_ferrari2014:LBD_chia2021:MS_imsgc2019:PD_blauwendraat2019:SCZ_pardinas2018:UC_delange2017
