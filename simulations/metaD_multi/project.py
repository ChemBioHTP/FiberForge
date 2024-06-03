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
    # Utility functions
######################
amino_acid_atoms = {
    'A': {'C': 3, 'H': 7, 'N': 1, 'O': 2},
    'R': {'C': 6, 'H': 14, 'N': 4, 'O': 2},
    'N': {'C': 4, 'H': 8, 'N': 2, 'O': 3},
    'D': {'C': 4, 'H': 7, 'N': 1, 'O': 4},
    'C': {'C': 3, 'H': 7, 'N': 1, 'O': 2},
    'Q': {'C': 5, 'H': 10, 'N': 2, 'O': 3},
    'E': {'C': 5, 'H': 9, 'N': 1, 'O': 4},
    'G': {'C': 2, 'H': 5, 'N': 1, 'O': 2},
    'H': {'C': 6, 'H': 9, 'N': 3, 'O': 2},
    'I': {'C': 6, 'H': 13, 'N': 1, 'O': 2},
    'L': {'C': 6, 'H': 13, 'N': 1, 'O': 2},
    'K': {'C': 6, 'H': 14, 'N': 2, 'O': 2},
    'M': {'C': 5, 'H': 11, 'N': 1, 'O': 2},
    'F': {'C': 9, 'H': 11, 'N': 1, 'O': 2},
    'P': {'C': 5, 'H': 9, 'N': 1, 'O': 1},
    'S': {'C': 3, 'H': 7, 'N': 1, 'O': 3},
    'T': {'C': 4, 'H': 9, 'N': 1, 'O': 3},
    'W': {'C': 11, 'H': 12, 'N': 2, 'O': 2},
    'Y': {'C': 9, 'H': 11, 'N': 1, 'O': 3},
    'V': {'C': 5, 'H': 11, 'N': 1, 'O': 2},
}

def get_n_atoms(sequence):
    return sum(sum(amino_acid_atoms[aa].values()) for aa in sequence)

######################
### Calculation Labels
######################
@Project.label
def created_dirs(job):
    return  (os.path.exists(job.path + '/0_build') and \
        os.path.exists(job.path + '/1_water_min') and \
        os.path.exists(job.path + '/2_water_eq') and \
        os.path.exists(job.path + '/3_sys_min') and \
        os.path.exists(job.path + '/4_sys_eq') and \
        os.path.exists(job.path + '/5_prod')
    )
@Project.label
def built_proteins(job):
    return job.isfile('0_build/seq.pdb')
@Project.label
def packed_proteins(job):
    return job.isfile('0_build/system.pdb')
@Project.label
def built_system_input(job):
    return job.isfile('0_build/tleap_system.in')
@Project.label
def built_system(job):
    return job.isfile('system.prmtop')
@Project.label
def finished_simulations(job):
    return job.isfile('5_prod/mdcrd.nc')
@Project.label
def ran_analysis(job):
    pass

######################
### Operations
######################
@Project.post(built_proteins)
@Project.post(packed_proteins)
@Project.post(built_system_input)
@Project.operation
def create_job_files(job):
    os.mkdir(job.path + '/0_build')
    os.mkdir(job.path + '/1_water_min')
    os.mkdir(job.path + '/2_water_eq')
    os.mkdir(job.path + '/3_sys_min')
    os.mkdir(job.path + '/4_sys_eq')
    os.mkdir(job.path + '/5_prod')
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
            '#load plumed',
            'module load GCC/8.2.0',
            'module load OpenMPI/3.1.4',
            'module load GROMACS/2019.2-PLUMED-2.5.2',
            "",
            "cp ../../in_files/min_wat.in 1_water_min/",
            "cp ../../in_files/md_wat.in 2_water_eq",
            "cp ../../in_files/min_sys.in 3_sys_min/",
            "cp ../../in_files/prod.in 5_prod/",
            "",
            "cd 1_water_min",
            "pmemd.cuda -O -i min_wat.in -o min_wat.out -p ../*.prmtop -c ../*.inpcrd -r min_wat.rst -ref ../*.inpcrd",
            "cp min_wat.rst ../2_water_eq",
            "cd ..",
            'echo "Finished step 1 "',
            "",
            "cd 2_water_eq",
            "pmemd.cuda -O -i md_wat.in -o md_wat.out -p ../*.prmtop -c ../1_water_min/*.rst -r md_wat.rst -ref ../1_water_min/*.rst",
            "cp md_wat.rst ../3_sys_min",
            "cd ..",
            "echo 'Finished step 2'",
            "",
            "cd 3_sys_min",
            "pmemd.cuda -O -i min_sys.in -o sys_min.out -p ../*.prmtop -c ../3_sys_min/*.rst -r sys_min.rst -ref ../3_sys_min/*.rst",
            "cp sys_min.rst ../4_sys_eq",
            "cd ..",
            "echo 'Finished step 3'",
            "",
            "cd 4_sys_eq",
            "pmemd.cuda -O -i heat.in -o heat.out -p ../*.prmtop -c ../4_sys_eq/*.rst -r heat.rst -ref ../4_sys_eq/*.rst",
            "cp heat.rst ../5_prod",
            "cd ..",
            "echo 'Finished step 4'",
            "",
            "cd 5_prod",
            "pmemd.cuda -O -i prod.in -o prod.out -p ../*.prmtop -c ../5_prod/heat.rst -r heat.rst -x mdcrd.nc",
            "echo 'Finished simulation'",
        ]
        for line in data:
            f.write(line + '\n')

    with open(job.path + '/4_sys_eq/heat.in', 'w') as f:
        data = [
            ' Heating System',
            '&cntrl',
            '   imin=0, nmropt=1,',
            '   ntx=1, irest=0,',
            '   ntpr=10000, ntwr=50000, ntwx=50000, iwrap=1,',
            '   ntf=2, ntb=1, cut=10.0, nsnb=20,',
            '   igb=0,',
            '   ibelly=0, ntr=1,',
            '   nstlim=250000, nscm=500, dt=0.002,',
            '   ntt=1, temp0=0.0, tempi=0.0, tautp=0.5',
            f"   ntc=2,restraintmask=':1-{len(job.sp.sequence)}',",
            '   restraint_wt=10.0,',
            '&end',
            '',
            "&wt type='REST', istep1=0, istep2=0, value1=1.0, value2=1.0, &end",
            "&wt type='TEMP0', istep1=0, istep2=250000, value1=0.0, value2=300, &end",
            "&wt type='END' &end",
        ]
        for line in data:
            f.write(line + '\n')
    
    atoms_per_molecule = get_n_atoms(job.sp.sequence)
    with open(job.path + '/5_prod/plumed.dat', 'w') as f:
        data = [
            '# define the center of the amyloid as the starting point',
            'm3: MOLECULES ...',
        ]

        data.append(f'  MOL1=1,{atoms_per_molecule}')
        for m in range(1, job.sp.n_molecules):
            data.append(f'   MOL{m+1}={m*atoms_per_molecule+1},{(m+1)*atoms_per_molecule}')
        data.append('...')
        data.append('')
        data.append('# define the distance between the centers of mass of the amyloids')
        data.extend([
            "s2: SMAC ...",
            "SPECIES=m3 LOWMEM",
            "KERNEL1={GAUSSIAN CENTER=0 SIGMA=0.480}  KERNEL2={GAUSSIAN CENTER=pi SIGMA=0.480}",
            'SWITCH={RATIONAL R_0=0.6}  MORE_THAN={RATIONAL R_0=0.7}  SWITCH_COORD={EXP R_0=4}',
            '...'
        ])
        data.append('')

        data.append('PRINT ARG=s2.* FILE=colvar STRIDE=10')

        for line in data:
            f.write(line + '\n')

@Project.pre(created_dirs)
@Project.post(built_proteins)
@Project.post(packed_proteins)
@Project.post(built_system_input)
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
@Project.pre(built_proteins)
@Project.post(packed_proteins)
@Project.post(built_system_input)
@Project.operation
def build_packmol_input(job):
    with open(job.path + '/0_build/packmol.inp', 'w') as f:
        data = [
            'tolerance 5.0',
            'filetype pdb',
            'output system.pdb',
            'structure seq_noH.pdb',
            f'   number {job.sp.n_molecules}',
            '   inside box 0. 0. 0. 150. 150. 150.',
            'end structure'
        ]
        for line in data:
            f.write(line + '\n')

    
    
@Project.pre(created_dirs)
@Project.pre(built_proteins)
@Project.pre(packed_proteins)
@Project.post(built_system_input)
@Project.operation
def build_tleap_system_input(job):
    # get charge of protein sequence
    total_charge = sum(aa_charge[aa] for aa in job.sp.sequence)
    if total_charge == 0:
        charge_line = ''
    else:
        charge_line = f'addIons agg Cl- {job.sp.n_molecules * total_charge}' if total_charge > 0 else f'addIons agg Na+ {job.sp.n_molecules * abs(total_charge)}'
    with open(job.path + '/0_build/tleap_system.in', 'w') as f:
        data = [
            'source leaprc.protein.ff14SB', 
            'source leaprc.gaff', 
            'source leaprc.water.tip3p',
            'agg = loadpdb system.pdb',
            'solvatebox agg TIP3PBOX 2',
            charge_line,
            'savepdb agg solvated_system.pdb',
            'saveamberparm agg system.prmtop system.inpcrd',
            'quit',
        ]
        for line in data:
            f.write(line + '\n')

@Project.pre(created_dirs)
@Project.post(packed_proteins)
@Project.post(built_system_input)
@Project.operation(cmd=True)
def run_protein(job):
    job.open()
    return f". ~/load_qz.sh; cd 0_build; tleap -f tleap_protein.in; pdb4amber --nohyd -i seq.pdb -o seq_noH.pdb;"
@Project.pre(created_dirs)
@Project.pre(built_proteins)
@Project.post(packed_proteins)
@Project.post(built_system_input)
@Project.operation(cmd=True)
def run_pack(job):
    job.open()
    return". ~/load_qz.sh; cd 0_build; packmol < packmol.inp"

@Project.pre(created_dirs)
@Project.pre(built_proteins)
@Project.pre(packed_proteins)
@Project.pre(built_system_input)
@Project.post(built_system)
@Project.operation(cmd=True)
def run_system_build(job):
    job.open()
    
    # packmol has an error that it doesn't add a TER at the end of the two protein sequences so we have to correct the pdb file before running it
    for i in range(job.sp.n_molecules):
        line_num = search_pattern_in_consecutive_lines(
            '0_build/system.pdb',
            pattern1=f'  {len(job.sp.sequence)}   ',
            pattern2=f'  1   ',
        )
        add_line_at_line_number(
            '0_build/system.pdb',
            line_number=line_num,
            line_content='TER'
        )

    return f". ~/load_qz.sh; cd 0_build; tleap -f tleap_system.in; cp system.prmtop system.inpcrd ..; cp solvated_system.pdb ../5_prod/"

# @Project.pre(created_dirs)
# @Project.pre(built_proteins)
# @Project.pre(packed_proteins)
# @Project.pre(built_system_input)
# @Project.pre(built_system)
# @Project.post(finished_simulations)
# @Project.operation(cmd=True)
# def run_simulation(job):
#     job.open()
#     return "sbatch run_simulation.sh"

if __name__ == "__main__":
    Project().main()
