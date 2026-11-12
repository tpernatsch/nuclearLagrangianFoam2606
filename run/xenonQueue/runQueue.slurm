#!/bin/bash
#SBATCH -n 32
#SBATCH -A eu
#SBATCH -J xenonQueue
#SBATCH -o xenonQueue.%j.out
# Runs a queue of xenon cases in sequence (runQueue.py). Submit from this folder:
#   sbatch [-n N] [-t HH:MM:SS] runQueue.slurm queueNextSteps.py [--only PATTERN ...]
# The OpenFOAM and repository environment come from the submitting shell (sbatch exports it):
#   source <openfoam>/etc/bashrc; source <repo>/etc/sourceme.sh
# Submitting again after a time-limit stop resumes the interrupted run and skips the finished ones.

cd "$SLURM_SUBMIT_DIR" || exit
export PATH=$FOAM_USER_APPBIN:$PATH

python3 runQueue.py "$@" --np "$SLURM_NTASKS"
