#!/bin/bash
#SBATCH --job-name=cg_nvt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=96:00:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --output=cg_nvt_%J.txt

WORKING_DIR=/home/nehilpkd/projects/biomatsims/cg_long/4_sys_eq
cd $WORKING_DIR

module load GCC/8.2.0
module load OpenMPI/3.1.4
module load GROMACS/2020

gmx mdrun -deffnm nvt -v

echo "The End"
