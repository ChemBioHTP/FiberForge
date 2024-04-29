#!/bin/bash
#SBATCH --job-name=meta
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=maxwell
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --output=%J.txt


source /home/shaoq1/bin/amber_env/amber-accre.sh


echo "here"
cd 1_water_min  
pmemd.cuda -O -i min_wat.in -o min_wat.out -p ../*.prmtop -c ../*.inpcrd -r min_wat.rst -ref ../*.inpcrd
echo "here"
cp min_wat.rst ../2_water_eq
cd ..
echo "Finished step 1 "

cd 2_water_eq  
pmemd.cuda -O -i md_wat.in -o md_wat.out -p ../*.prmtop -c ../1_water_min/*.rst -r md_wat.rst -ref ../1_water_min/*.rst
cp md_wat.rst ../3_sys_min
cd ..
echo "Finished step 2"

cd 3_sys_min  
pmemd.cuda -O -i min_sys.in -o sys_min.out -p ../*.prmtop -c ../3_sys_min/*.rst -r sys_min.rst -ref ../3_sys_min/*.rst
cp sys_min.rst ../4_sys_eq
cd ..
echo "Finished step 3"

cd 4_sys_eq  
pmemd.cuda -O -i heat.in -o heat.out -p ../*.prmtop -c ../4_sys_eq/*.rst -r heat.rst -ref ../4_sys_eq/*.rst
cp heat.rst ../5_prod
cd ..
echo "Finished step 4"

cd 5_prod
pmemd.cuda -O -i prod.in -o prod.out -p ../*.prmtop -c ../5_prod/heat.rst -r heat.rst -x mdcrd.nc
echo "The End"
