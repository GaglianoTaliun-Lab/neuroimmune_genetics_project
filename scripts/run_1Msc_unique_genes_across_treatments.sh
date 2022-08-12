#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=4:00:00
#SBATCH --job-name=get_unique_genes_across_treatments
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=5G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

Rscript 1Msc_unique_genes_across_treatments.R
