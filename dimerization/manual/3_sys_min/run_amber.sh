#!/bin/bash
#SBATCH --job-name=foo_amb
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=maxwell
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=05:00:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --output=test_%J.txt


source /home/shaoq1/bin/amber_env/amber-accre.sh

pmemd.cuda -O -i min_sys.in -o sys_min.out -p ../*.prmtop -c ../3_sys_min/*.rst -r sys_min.rst -ref ../3_sys_min/*.rst

echo "The End"
