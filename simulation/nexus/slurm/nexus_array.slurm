#!/bin/bash
# One NEXUS batch per array task.
#   sbatch --array=0-(N-1) --output=OUTDIR/logs/%x.%A_%a.out --error=OUTDIR/logs/%x.%A_%a.err \
#          nexus_array.slurm SAMPLE EVENTS OUTDIR [SEED_OFFSET]
# SAMPLE: bb0nu | electron | bi214. Batch number = array task; event numbers start at
# batch*EVENTS; see run_nexus_batch.sh.
#SBATCH --job-name=nexus_itaca
#SBATCH --partition=general
#SBATCH --qos=serial
#SBATCH --constraint=icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
/scratch/jjgomez/nexus/slurm/run_nexus_batch.sh "$1" "$SLURM_ARRAY_TASK_ID" "$2" "$3" "${4:-0}"
