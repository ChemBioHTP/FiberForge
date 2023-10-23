#!/bin/bash
#SBATCH --job-name=foo_amb
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=pascal
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=20:00:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --output=test_%J.txt


source /home/shaoq1/bin/amber_env/amber-accre.sh

#load plumed
module load GCC/8.2.0
module load OpenMPI/3.1.4
module load GROMACS/2019.2-PLUMED-2.5.2

pmemd.cuda -O -i prod.in -o prod.out -p ../*.prmtop -c ../5_meta_COM/heat.rst -r heat.rst


echo "The End"
