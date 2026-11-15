#!/bin/bash
#SBATCH --job-name=nexus_itaca_test
#SBATCH --partition=general
#SBATCH --qos=test
#SBATCH --constraint=icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=00:10:00
#SBATCH --array=0-1
#SBATCH --output=/scratch/jjgomez/NAUSICAA-runs/nexus_test/%x.%A_%a.out
#SBATCH --error=/scratch/jjgomez/NAUSICAA-runs/nexus_test/%x.%A_%a.err
SAMPLES=(bb0nu electron)
/scratch/jjgomez/nexus/slurm/run_nexus_batch.sh ${SAMPLES[$SLURM_ARRAY_TASK_ID]} 0 500 /scratch/jjgomez/NAUSICAA-runs/nexus_test
