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
    ### Utility functions
######################
def find_repeating_amino_acids(file_path):
    import re

    def remove_numbers(input_string):
        return re.sub(r'\d+', '', input_string)

    with open(file_path, 'r') as file:
        lines = file.readlines()
        file.close()
        
    amino_acid_code = remove_numbers(lines[2].split()[0])
    start = 1
    end = 1
    repeats = []
    for line_number, line in enumerate(lines[3:], start=1): # start at line 3
        new_amino_acid_code = remove_numbers(line.split()[0])
        if new_amino_acid_code != amino_acid_code:
            amino_acid_code = new_amino_acid_code
            repeats.append((start, end))
            start = line_number
        end = line_number
    return repeats
        

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
def finished_nvt_simulation(job):
    return job.isfile('3_sys_eq/nvt.gro')
@Project.label
def finished_simulations(job):
    return job.isfile('4_prod/prod.gro')
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
            "cp ../../in_files/npt.mdp 4_prod",
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
            "",
            "cd 4_prod",
            "gmx grompp -f npt.mdp -c ../3_sys_eq/nvt.gro -p ../topol.top -o prod.tpr -maxwarn 1",
            "gmx mdrun -plumed plumed.dat -deffnm prod",
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

@Project.pre(built_proteins)
@Project.pre(finished_nvt_simulation)
@Project.post(finished_simulations)
@Project.operation
def create_plumed(job):
    job.open()
    repeats = find_repeating_amino_acids(job.fn('3_sys_eq/nvt.gro'))
    data = []
    data.append('#treat each molecule as whole')
    whole = "WHOLEMOLECULES "
    whole += f"ENTITY0=1-{repeats[len(job.sp.sequence)-1][1]} "
    for i in range(1, job.sp.n_molecules):
        whole += f"ENTITY{i}={repeats[(i)*len(job.sp.sequence)][0]}-{repeats[(i+1)*len(job.sp.sequence)-1][1]} "
    data.append(whole)
    data.extend([
        '# define the center of the amyloid as the starting point',
        'm3: MOLECULES ...',
    ])
    data.append(f'  MOL1=1,{repeats[len(job.sp.sequence)-1][1]}')
    for m in range(2, job.sp.n_molecules):
        data.append(f'   MOL{m}={repeats[(m-1)*len(job.sp.sequence)][0]},{repeats[(m)*len(job.sp.sequence)-1][1]}')
    data.append('...')
    data.append('')
    data.append('# define the collective variables')
    data.extend([
        "s2: SMAC ...",
        "SPECIES=m3 LOWMEM",
        "KERNEL1={GAUSSIAN CENTER=0 SIGMA=0.480}  KERNEL2={GAUSSIAN CENTER=pi SIGMA=0.480}",
        'SWITCH={RATIONAL R_0=0.6}  MORE_THAN={RATIONAL R_0=0.7}  SWITCH_COORD={EXP R_0=4}',
        '...'
    ])
    data.append('')
    data.append('# define the metadynamics')
    data.append('DUMPMULTICOLVAR DATA=s2 FILE=MULTICOLVAR.xyz STRIDE=1000')

    with open(job.fn('4_prod/plumed.dat'), 'w') as f:
        for line in data:
            f.write(line + '\n')

@Project.pre(created_dirs)
@Project.pre(built_proteins)
@Project.post(finished_simulations)
@Project.post(finished_nvt_simulation)
@Project.operation(cmd=True)
def run_simulation(job):
    job.open()
    return "sbatch run_simulation.sh"

if __name__ == "__main__":
    Project().main()
