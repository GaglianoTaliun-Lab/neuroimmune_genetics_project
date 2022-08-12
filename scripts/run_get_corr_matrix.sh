#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=2:00:00
#SBATCH --job-name=get_sample_overlap
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

# arguments:
# arg 1 = lava run version (e.g. v1, v2, v3)
# arg 2-x = phenotypes to include in the sample overlap matrix

Rscript get_corr_matrix.R
