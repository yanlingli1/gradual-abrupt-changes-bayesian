#!/bin/bash -l 
#SBATCH --nodes=1 
#SBATCH --ntasks=2
#SBATCH --mem=72GB 
#SBATCH --time=48:00:00 
#SBATCH --account=open
#SBATCH --array=99  # use this to submit multiple jobs
#SBATCH --output=submit_rsgcm.out # comment this out if you want the output file for each job     

# make sure to use R in the conda environment
R_EXEC=/storage/work/xjx5093/.conda/envs/jagstest/bin/R

N=100
O=200
shift=2
r1=$SLURM_ARRAY_TASK_ID

module load anaconda 
#source /swst/apps/anaconda3/2021.05_gcc-8.5.0/etc/profile.d/conda.sh
conda activate jagstest

# make sure to use the conda environment
export PATH=/storage/work/xjx5093/.conda/envs/jagstest/bin:$PATH

# make sure to library rjags from the conda environment
export R_LIBRARY_LOAD_PATH=/storage/work/xjx5093/.conda/envs/jagstest/lib/R/library

cd /storage/work/xjx5093/GOHIAR_RS/RS-GCM_SparseTransition

eval "${R_EXEC} CMD BATCH --no-save --no-restore '--args N=${N} O=${O} shift=${shift} r1=${r1}' Bayesian_RSTVP_gcm_final.R Bayesian_RSTVP_gcm_final_shift${shift}_N${N}O${O}.Rout" 

conda deactivate

