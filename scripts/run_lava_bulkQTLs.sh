#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=1-00:00:00
#SBATCH --job-name=lava_sign_bulk_qtls_array
#SBATCH --output=slurm-%x-array.out
#SBATCH --error=slurm-%x-array.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=5G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

#################### LAVA with QTLs:
# Rscript lava_bulkQTL_target_loci_A.R
Rscript lava_bulkQTL_target_loci_significant.R
