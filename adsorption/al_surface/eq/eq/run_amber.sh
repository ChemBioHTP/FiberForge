#!/bin/bash
#SBATCH --job-name=foo_amb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=10:00:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --output=test_%J.txt

WORKING_DIR=/home/nehilpkd/projects/biomatsims/adsorption/al_surface/eq
cd $WORKING_DIR

source /home/shaoq1/bin/amber_env/amber-accre.sh

sander -O -i mfp_al_heat.in -o mfp_al_heat.out -p ../mfp_al.prmtop -c mfp_al_min.rst -r mfp_al_heat.rst -x heat.netcdf -ref mfp_al_min.rst


echo "The End"
