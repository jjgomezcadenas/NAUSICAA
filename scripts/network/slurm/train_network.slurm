#!/bin/bash
# One training of the classification network on one GPU. Hyperion cluster.
#   sbatch --job-name=NAME --qos=regular --time=24:00:00 --output=LOGDIR/%x.%j.out --error=LOGDIR/%x.%j.err \
#          train_network.slurm SIGNAL.h5 BACKGROUND.h5 OUTPUT_DIR --input=itrack [train_network.py options]
# Paths are relative to the repository. Any GPU type of the partition is accepted.
#SBATCH --job-name=itaca_network
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --gres=gpu:1
set +e
type module >/dev/null 2>&1 || source /etc/profile
set -eo pipefail
module purge
module load PyTorch/2.1.2-foss-2021a-CUDA-12.1.1
# The module has no h5py: the environment adds it on top of the module's Python (see requirements-network.txt)
PY=${NAUSICAA_NETWORK_PYTHON:-/scratch/jjgomez/network-venv/bin/python}
cd "${NAUSICAA_REPO:-/scratch/jjgomez/NAUSICAA}"
echo "node $(hostname) commit $(git rev-parse HEAD) gpu $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | head -1)"
$PY -c "import torch,h5py,numpy;print('torch',torch.__version__,'cuda',torch.cuda.is_available(),'h5py',h5py.__version__,'numpy',numpy.__version__)"
$PY scripts/network/train_network.py "$@"
