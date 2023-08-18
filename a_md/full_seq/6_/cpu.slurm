#!/bin/bash
#SBATCH --mail-user=kieran.d.nehil-puleo@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=1000M
#SBATCH --output=cpu_accel_md.out


WORKING_DIR=/home/nehilpkd/projects/foo/6_
cd $WORKING_DIR


source /home/shaoq1/bin/amber_env/amber-accre.sh

pmemd -O -i amd.in -o amd.out -p ../*.prmtop -c ../6_/eq.rst -r prod.nc
