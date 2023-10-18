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

WORKING_DIR=/home/nehilpkd/projects/biomatsims/a_md/mfp/3_
cd $WORKING_DIR

source /home/shaoq1/bin/amber_env/amber-accre.sh

pmemd.cuda -O -i heat.in -o heat.out -p ../*.prmtop -c ../3_/*.rst -r heat.rst -ref ../3_/*.rst

echo "The End"
