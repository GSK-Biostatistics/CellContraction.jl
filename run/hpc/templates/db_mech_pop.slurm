#!/bin/bash

#SBATCH --job-name=cellcontr
#SBATCH --partition=cpu
#SBATCH --array=${IDX}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${NWORKERS}
#SBATCH --time=01:00:00

${JULIA} --project=${PROJECT} ${SCRIPT}
