#!/bin/bash
#SBATCH -p test.q
#SBATCH -n 1
module purge
module load likwid/4.1
module load slurm/current
