#!/bin/bash
#SBATCH --job-name=dd03d5ae615caf735105b01b56166598
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=pascal
#SBATCH --gres=gpu:1
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mem=10G
#SBATCH --output=job.log
source /home/shaoq1/bin/amber_env/amber-accre.sh
module load GCC/8.2.0
module load OpenMPI/3.1.4
module load GROMACS/2020
source ~/.bashrc
conda activate biomatsims

cp ../../in_files/minim.mdp 1_min
cp ../../in_files/nvt.mdp 2_eq_nvt
cp ../../in_files/npt.mdp 3_eq_npt

cd 1_min
gmx grompp -f minim.mdp -c ../0_preprocess/solvated_ions.gro -p ../0_preprocess/topol.top -o em.tpr
gmx mdrun -v -deffnm em
cd ..
cd 2_eq_nvt
gmx grompp -f nvt.mdp -c ../1_min/em.gro -p ../0_preprocess/topol.top -r ../1_min/em.gro -o nvt.tpr
gmx mdrun -deffnm nvt
cd ..
cd 3_eq_npt
gmx grompp -f npt.mdp -c ../2_eq_nvt/nvt.gro -p ../0_preprocess/topol.top -r ../2_eq_nvt/nvt.gro -o npt.tpr
gmx mdrun -deffnm npt
cd ..
echo 'Finished equilibration'
