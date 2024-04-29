#!/bin/bash
#SBATCH --job-name=smd_710d5fb75ec5a28ce9a4d53c40ecd5c0
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=pascal
#SBATCH --gres=gpu:1
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mem=10G
#SBATCH --output=pull.log
source /home/shaoq1/bin/amber_env/amber-accre.sh
module load GCC/8.2.0
module load OpenMPI/3.1.4
module load GROMACS/2020
source ~/.bashrc
conda activate biomatsims

cd 4_smd
gmx mdrun -deffnm pull
echo 'Finished pulling'
