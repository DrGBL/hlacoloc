#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=25G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-999
#SBATCH --output=/path/to/array_test_%A_%a.out
#SBATCH --error=/path/to/array_test_%A_%a.error


cd /path/to/simulation/

anc_vector=( serology afr amr eas sas )
anc_var=${anc_vector[ $(( SLURM_ARRAY_TASK_ID / 200 )) ]}

iter_var=$(( (SLURM_ARRAY_TASK_ID ) % 200 +1 ))

#binary L=10
Rscript slurm_simulation_script.R \
  --anc ${anc_var} \
  --ids /path/to/ukb.${anc_var}IDsPCA.txt \
  --call_folder /path/to/calls/ \
  --af_folder /path/to/simulation/ \
  --bim_folder /path/to/simulation/ \
  --out /path/to/simulation/slurm_sims_binary/ \
  --n_min_alleles 10 \
  --pheno2_mean_shared_factor 0 \
  --pheno2_variance_shared_factor 1 \
  --gen_var 0.1 \
  --freq_threshold 0.01 \
  --iteration ${iter_var}_r2_10perc \
  --plot_susie FALSE \
  --continuous_trait FALSE \
  --n_clone 10

#quantitative L=10
Rscript slurm_simulation_script.R \
  --anc ${anc_var} \
  --ids /path/to/ukb.${anc_var}IDsPCA.txt \
  --call_folder /path/to/calls/ \
  --af_folder /path/to/simulation/ \
  --bim_folder /path/to/simulation/ \
  --out /path/to/simulation/slurm_sims/ \
  --n_min_alleles 10 \
  --pheno2_mean_shared_factor 0 \
  --pheno2_variance_shared_factor 1 \
  --gen_var 0.1 \
  --freq_threshold 0.01 \
  --iteration ${iter_var}_r2_10perc \
  --plot_susie FALSE \
  --continuous_trait TRUE \
  --n_clone 1

#binary L=20
Rscript slurm_simulation_script.R \
  --anc ${anc_var} \
  --ids /path/to/ukb.${anc_var}IDsPCA.txt \
  --call_folder /path/to/calls/ \
  --af_folder /path/to/simulation/ \
  --bim_folder /path/to/simulation/ \
  --out /path/to/simulation/slurm_sims_binary_L20/ \
  --n_min_alleles 10 \
  --pheno2_mean_shared_factor 0 \
  --pheno2_variance_shared_factor 1 \
  --gen_var 0.1 \
  --freq_threshold 0.01 \
  --iteration ${iter_var}_r2_10perc \
  --plot_susie FALSE \
  --continuous_trait FALSE \
  --n_clone 10 \
  --susie_L 20

#quantitative L=20
Rscript slurm_simulation_script.R \
  --anc ${anc_var} \
  --ids /path/to/ukb.${anc_var}IDsPCA.txt \
  --call_folder /path/to/calls/ \
  --af_folder /path/to/simulation/ \
  --bim_folder /path/to/simulation/ \
  --out /path/to/simulation/slurm_sims_L20/ \
  --n_min_alleles 10 \
  --pheno2_mean_shared_factor 0 \
  --pheno2_variance_shared_factor 1 \
  --gen_var 0.1 \
  --freq_threshold 0.01 \
  --iteration ${iter_var}_r2_10perc \
  --plot_susie FALSE \
  --continuous_trait TRUE \
  --n_clone 1 \
  --susie_L 20