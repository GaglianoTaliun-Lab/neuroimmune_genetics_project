#!/bin/bash
#SBATCH --account=def-gsarah
#SBATCH --time=1:00:00
#SBATCH --array=1-9
#SBATCH --job-name=munge_sumstats
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3G

module load StdEnv/2020
module load python/2.7
module load scipy-stack/2020a

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --upgrade pip
# pip install bitarray
# pip install pybedtools==0.7.10
pip install contextlib2==0.2
pip install importlib_metadata==1.5.0
pip install jsonschema==2.6.0
pip install six==1.16.0+computecanada
pip install bitarray==0.8.0

project_dir="/home/fridald4/projects/def-gsarah/fridald4/neuroimmune_genetics_project"

pip install -r ${project_dir}/tools/ldsc/CC_requirements.txt

study=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/GWAS_summary_statistics/list_of_GWAS_studies.txt | awk '{print $1}')
total_N=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/GWAS_summary_statistics/list_of_GWAS_studies.txt | awk '{print $2}')

study_noprox=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/GWAS_summary_statistics/list_of_GWAS_studies_no_proxies.txt | awk '{print $1}')
total_N_noprox=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/GWAS_summary_statistics/list_of_GWAS_studies_no_proxies.txt | awk '{print $2}')

python ${project_dir}/tools/ldsc/munge_sumstats.py \
	--sumstats ${project_dir}/GWAS_summary_statistics/tsv_for_munge_stats/${study}.tsv \
	--out ${project_dir}/GWAS_summary_statistics/ldsc_sumstats/${study} \
	--N ${total_N} \
	--chunksize 500000 \ 
	--merge-alleles ${project_dir}/reference_data/w_hm3.snplist

python ${project_dir}/tools/ldsc/munge_sumstats.py \
        --sumstats ${project_dir}/GWAS_summary_statistics/tsv_for_munge_stats/${study_noprox}.tsv \
        --out ${project_dir}/GWAS_summary_statistics/ldsc_sumstats/${study_noprox} \
        --N ${total_N_noprox} \
        --chunksize 500000 \ 
        --merge-alleles ${project_dir}/reference_data/w_hm3.snplist
