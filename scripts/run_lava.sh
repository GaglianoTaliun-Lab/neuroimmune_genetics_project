#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=1-0:00:00
#SBATCH --job-name=lava_no_proxies
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

##################### LAVA for GWAS traits

# Rscript lava.R
Rscript lava_no_proxies.R
