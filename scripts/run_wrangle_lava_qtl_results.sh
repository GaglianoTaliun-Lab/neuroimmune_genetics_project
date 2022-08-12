#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=0:30:00
#SBATCH --job-name=wrangle-lava-qtl
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

# Rscript wrangle_lava_qtl_results.R
# Rscript wrangle_lava_bulkQTL_results.R
# Rscript wrangle_lava_onek1k_results.R
Rscript wrangle_lava_pqtls_results.R
# Rscript wrangle_lava_1MscBloodNL_treatments_results.R
