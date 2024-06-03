import os
import pathlib
from copy import deepcopy

import numpy as np
from flow import FlowProject, staticlabel, aggregator
from flow.environment import DefaultSlurmEnvironment
import signac
from textwrap import wrap
import signac
from fiberForge.utils import *

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)
warnings.simplefilter(action='ignore')


class Project(FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""
    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()

class Accre(DefaultSlurmEnvironment):
    """Subclass of DefaultPBSEnvironment for VU's Rahman cluster."""

    # template = "rahman.sh"

######################
### Calculation Labels
######################
@Project.label
def created_dirs(job):
    return  (os.path.exists(job.path + '/0_build') and \
        os.path.exists(job.path + '/1_ions') and \
        os.path.exists(job.path + '/2_sys_min') and \
        os.path.exists(job.path + '/3_sys_eq') and \
        os.path.exists(job.path + '/4_prod')
    )
@Project.label
def built_proteins(job):
    return job.isfile('0_build/seq.pdb')
@Project.label
def finished_simulations(job):
    return job.isfile('3_sys_eq/nvt.xtc')
@Project.label
def ran_analysis(job):
    pass

######################
### Operations
######################
@Project.post(built_proteins)
@Project.operation
def create_job_files(job):
    os.mkdir(job.path + '/0_build')
    # os.mkdir(job.path + '/1_water_eq')
    os.mkdir(job.path + '/1_ions')
    os.mkdir(job.path + '/2_sys_min')
    os.mkdir(job.path + '/3_sys_eq')
    os.mkdir(job.path + '/4_prod')
    with open(job.path + '/run_simulation.sh', 'w') as f:
        data = [
            "#!/bin/bash",
            f"#SBATCH --job-name={job.id}",
            "#SBATCH --account=csb_gpu_acc",
            "#SBATCH --partition=pascal",  
            "#SBATCH --gres=gpu:1",
            "#SBATCH --no-requeue",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks-per-node=1",
            "#SBATCH --time=72:00:00",
            "#SBATCH --mem=10G",            
            "#SBATCH --output=job.log",
            ""
            "source /home/shaoq1/bin/amber_env/amber-accre.sh",
            "module load GCC/8.2.0",
            "module load OpenMPI/3.1.4",
            "module load GROMACS/2019.2-PLUMED-2.5.2",
            "source ~/.bashrc",
            "conda activate fiberForge",
            "",
            "cp ../../in_files/12nm_water_eq.gro 0_build",
            "cp ../../in_files/martinize.py 0_build",
            "cp ../../in_files/martini_v2.0_ions.itp 0_build",
            "cp ../../in_files/martini_v2.2.itp 0_build",
            "cp ../../in_files/ions.mdp 1_ions",
            "cp ../../in_files/water_min.mdp 2_sys_min",
            "cp ../../in_files/nvt.mdp 3_sys_eq",
            "",
            "cd 0_build",
            f"python martinize.py -f seq.pdb -o topol.top -x seq_cg.pdb -ff martini22",
            f"gmx insert-molecules -box 12 12 12 -nmol {job.sp.n_molecules} -ci seq_cg.pdb -radius 0.4 -o seq_{job.sp.n_molecules}_cg.gro -try 1000",
            f"gmx solvate -cp seq_{job.sp.n_molecules}_cg.gro -cs 12nm_water_eq.gro -radius 0.21 -o seq_{job.sp.n_molecules}_cg_water.gro -p topol.top",
            "sed -i '1d' topol.top",
            "sed -i '1s/^/#include \"martini_v2.0_ions.itp\"/' topol.top",
            "sed -i '1s/^/\\n/' topol.top",
            "sed -i '1s/^/#include \"martini_v2.2.itp\"/' topol.top",
            "sed -i 's/\(Protein[[:space:]]*\)\([0-9A-Za-z]\)[[:space:]]*\([0-9]*\)/\\1\\2\\n\\3/' topol.top",
            f"sed -i 's/Protein[[:space:]]*1$/Protein          {job.sp.n_molecules}/' topol.top",
            "cp topol.top Protein.itp martini_v2.2.itp martini_v2.0_ions.itp ..",
            "cd ..",
            'echo "Finished step 0 "',
            "",
            "cd 1_ions",
            f"gmx grompp -f ions.mdp -c ../0_build/seq_{job.sp.n_molecules}_cg_water.gro -p ../topol.top -o ions.tpr",
            f"echo W | gmx genion -s ions.tpr -o seq_{job.sp.n_molecules}_cg_water_ions.gro -p ../topol.top -pname NA+ -nname CL- -neutral | 13",
            "cd ..",
            'echo "Finished step 1 "',
            "",
            "cd 2_sys_min",
            f"gmx grompp -f water_min.mdp -c ../1_ions/seq_{job.sp.n_molecules}_cg_water_ions.gro -p ../topol.top -o em.tpr -maxwarn 1",
            "gmx mdrun -deffnm em",
            "cd ..",
            "echo 'Finished step 2'",
            "",
            "cd 3_sys_eq",
            "gmx grompp -f nvt.mdp -c ../2_sys_min/em.gro -p ../topol.top -o nvt.tpr -maxwarn 1",
            "gmx mdrun -deffnm nvt",
            "cd ..",
            "echo 'Finished step 3'",
            # "",
            # "cd 5_prod",
            # "pmemd.cuda -O -i prod.in -o prod.out -p ../*.prmtop -c ../5_prod/heat.rst -r heat.rst -x mdcrd.nc",
            # "echo 'Finished simulation'",
        ]
        for line in data:
            f.write(line + '\n')

@Project.pre(created_dirs)
@Project.post(built_proteins)
@Project.operation
def build_tleap_input(job):
    with open(job.path + '/0_build/tleap_protein.in', 'w') as f:
        data = [
            'source leaprc.protein.ff14SBonlysc', 
            f'mol = sequence {{{single_to_triple(job.sp.sequence)}}}',
            'savepdb mol seq.pdb',
            'quit',
        ]
        for line in data:
            f.write(line + '\n')
    

@Project.pre(created_dirs)
@Project.post(built_proteins)
@Project.operation(cmd=True)
def run_protein(job):
    job.open()
    return f". ~/load_qz.sh; cd 0_build; tleap -f tleap_protein.in;"

@Project.pre(created_dirs)
@Project.pre(built_proteins)
@Project.post(finished_simulations)
@Project.operation(cmd=True)
def run_simulation(job):
    job.open()
    return "sbatch run_simulation.sh"

if __name__ == "__main__":
    Project().main()
