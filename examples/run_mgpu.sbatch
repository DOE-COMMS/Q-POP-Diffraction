#!/bin/bash
# Multi-GPU diffraction of one FerroDyn snapshot with Slurm.
#
# Submit from the FerroDyn run directory, next to input.json and parameter.atom.in.
# Adapt the account, partition, GPU options and environment to your system, and use as
# many GPUs as `qpop-diffraction --mgpu` asks for (e.g. 8 GPUs of 40 GB for 4000x4000x68).
#
#SBATCH -A <account>
#SBATCH -p <gpu partition>
#SBATCH -N 2
#SBATCH -t 00:30:00
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH -o diffraction_%j.out

set -euo pipefail

QPOP_DIFFRACTION=${QPOP_DIFFRACTION:-$HOME/Q-POP-Diffraction}

# The compiler, CUDA, MPI, HDF5 and HeFFTe used for the multi-GPU build, and GPU-aware MPI
# (for Cray MPICH: MPICH_GPU_SUPPORT_ENABLED=1)
# module load ...
export MPICH_GPU_SUPPORT_ENABLED=1

# Setup in Python: writes <output_file>_params.h5 and <output_file>.h5 without /intensity
qpop-diffraction -i input.json -s system_setting.json -a parameter.atom.in --mgpu

# One MPI rank per GPU; each rank selects its GPU itself
OUTPUT=$(python3 -c "import json; print(json.load(open('input.json'))['output_file'])")
srun --ntasks-per-node=4 --gpus-per-node=4 --gpu-bind=none \
    "$QPOP_DIFFRACTION/bin/qpop_diffraction_mgpu" "${OUTPUT}_params.h5"
