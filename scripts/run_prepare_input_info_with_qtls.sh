#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=1:00:00
#SBATCH --job-name=prepare_input_info_with_qtls
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

# Rscript prepare_input_info_with_qtls.R
# Rscript prepare_input_info_with_onek1k.R
# Rscript prepare_input_info_with_bulkQTLs.R
# Rscript prepare_input_info_1Msc_treatments.R
Rscript prepare_input_info_pqtls.R
