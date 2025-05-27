#!/bin/bash
#SBATCH -A e30514
#SBATCH -p short
#SBATCH --job-name="diffusion"
#SBATCH -n 4
#SBATCH -t 00:01:30
#SBATCH --mem=1G
#SBATCH --mail-user=piercedillon2026@u.northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL

cd $SLURM_SUBMIT_DIR
module load mpi/openmpi-4.1.7-gcc-10.4.0
mpirun -np 4 diffusion 512 1.95 1e-9 30000
