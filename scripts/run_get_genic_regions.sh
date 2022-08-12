#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=2:00:00
#SBATCH --job-name=get_genic_regions_GWS
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

# Rscript get_genic_regions_GWS.R
# Rscript get_genic_onek1k_regions_GWS.R
# Rscript get_genic_1Msc_treatments_regions.R
Rscript get_genic_pqtls_regions.R
sed -i 's/"//g' pqtls_filtered.loci # only run it for the pQTLs, since I cant seem to remove the quoting marks on some genes.
