#!/bin/bash
#SBATCH --job-name=foo_amb
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=pascal
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=05:00:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --output=test_%J.txt

WORKING_DIR=/home/nehilpkd/projects/biomatsims/aggregation/0_
cd $WORKING_DIR


source /home/shaoq1/bin/amber_env/amber-accre.sh

pmemd.cuda -O -i min_wat.in -o min_wat.out -p ../*.prmtop -c ../*.inpcrd -r min_wat.rst -ref ../*.inpcrd

echo "The End"
