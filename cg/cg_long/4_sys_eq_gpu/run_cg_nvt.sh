#!/bin/bash
#SBATCH --job-name=nvt_gpu
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=maxwell
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --output=nvt_%J.txt

WORKING_DIR=/home/nehilpkd/projects/biomatsims/cg_long/4_sys_eq_gpu
cd $WORKING_DIR

module load GCC/8.2.0
module load OpenMPI/3.1.4
module load GROMACS/2020

gmx mdrun -deffnm nvt -v

echo "The End"
