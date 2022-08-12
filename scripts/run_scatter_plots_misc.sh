#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=1:00:00
#SBATCH --job-name=scatter-plots
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=3G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

Rscript scatter_plots_misc.R
