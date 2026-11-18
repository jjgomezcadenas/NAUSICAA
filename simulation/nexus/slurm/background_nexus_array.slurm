#!/bin/bash
# One NEXUS batch of a radiogenic background run (decays in the copper shield) per array task.
#   sbatch --array=0-(N-1) --output=OUTDIR/logs/%x.%A_%a.out --error=OUTDIR/logs/%x.%A_%a.err \
#          background_nexus_array.slurm NAME ISOTOPE EVENTS OUTDIR [SEED_OFFSET]
# NAME: folder of macros/background (II_L1500_ICS20, IE_L1500_ICS20), written from
# the shielding cards by scripts/background/write_background_macros.jl.
# ISOTOPE: bi214 | tl208. Batch number = array task; event numbers start at
# batch*EVENTS; see run_iitaca_nexus_batch.sh (MACRO_SET=background).
#SBATCH --job-name=nexus_background
#SBATCH --partition=general
#SBATCH --exclude=hyperion-219   # this node fails to resolve the user: tasks there are requeued and held ("user env retrieval failed", 2026-11-17)
#SBATCH --qos=serial
#SBATCH --constraint=icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=04:00:00
MACRO_SET=background /scratch/jjgomez/nexus/slurm/run_iitaca_nexus_batch.sh "$1" "$2" "$SLURM_ARRAY_TASK_ID" "$3" "$4" "${5:-0}"
