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
        os.path.exists(job.path + '/1_water_min') and \
        os.path.exists(job.path + '/2_water_eq') and \
        os.path.exists(job.path + '/3_sys_min') and \
        os.path.exists(job.path + '/4_sys_eq') and \
        os.path.exists(job.path + '/5_meta')
    )
@Project.label
def built_proteins(job):
    return job.isfile('0_build/seq1.pdb') and job.isfile('0_build/seq2.pdb')
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
    return job.isfile('5_meta/mdcrd.nc')
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
    os.mkdir(job.path + '/5_meta')
    with open(job.path + '/run_simulation.sh', 'w') as f:
        data = [
            "#!/bin/bash",
            f"#SBATCH --job-name=meta_{job.id}",
            "#SBATCH --account=csb_gpu_acc",
            "#SBATCH --partition=pascal",  
            "#SBATCH --gres=gpu:1",
            "#SBATCH --no-requeue",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks-per-node=1",
            "#SBATCH --time=48:00:00",
            "#SBATCH --mem=10G",            
            "#SBATCH --output=job.log",
            ""
            "source /home/shaoq1/bin/amber_env/amber-accre.sh",
            "module load GCC/8.2.0",
            "module load OpenMPI/3.1.4",
            "module load GROMACS/2019.2-PLUMED-2.5.2",
            "",
            "cp ../../in_files/min_wat.in 1_water_min/",
            "cp ../../in_files/md_wat.in 2_water_eq",
            "cp ../../in_files/min_sys.in 3_sys_min/",
            "cp ../../in_files/prod.in 5_meta/",
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
            "cp heat.rst ../5_meta",
            "cd ..",
            "echo 'Finished step 4'",
            "",
            "cd 5_meta",
            "pmemd.cuda -O -i prod.in -o prod.out -p ../*.prmtop -c ../5_meta/heat.rst -r heat.rst -x mdcrd.nc",
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
            '   ntpr=10000, ntwr=500, ntwx=500, iwrap=1,',
            '   ntf=2, ntb=1, cut=10.0, nsnb=20,',
            '   igb=0,',
            '   ibelly=0, ntr=1,',
            '   nstlim=250000, nscm=500, dt=0.002,',
            '   ntt=1, temp0=0.0, tempi=0.0, tautp=0.5',
            f"   ntc=2,restraintmask=':1-{len(job.sp.amyloids[0])+len(job.sp.linkers[0]) + len(job.sp.amyloids[1])+len(job.sp.linkers[1])}',",
            '   restraint_wt=10.0,',
            '&end',
            '',
            "&wt type='REST', istep1=0, istep2=0, value1=1.0, value2=1.0, &end",
            "&wt type='TEMP0', istep1=0, istep2=250000, value1=0.0, value2=300, &end",
            "&wt type='END' &end",
        ]
        for line in data:
            f.write(line + '\n')
    with open(job.path + '/5_meta/plumed.dat', 'w') as f:
        data = [
            'MOLINFO STRUCTURE=solvated_system.pdb',
            '',
            '# define CV',
            f'ab: ANTIBETARMSD RESIDUES={job.sp.amyloid_location_index}-{len(job.sp.amyloids[0])},{job.sp.amyloid_location_index}-{len(job.sp.amyloids[1])} STRANDS_CUTOFF=1',
            '',
            '# metadynamics argument',
            'metad: METAD ARG=ab PACE=500 HEIGHT=0.2 SIGMA=0.05 BIASFACTOR=8 GRID_MIN=-0.1 GRID_MAX=3.0 GRID_BIN=500 FILE=HILLS TEMP=300',
            '',
            '# monitor distance',
            'PRINT STRIDE=500 ARG=ab,metad.bias FILE=COLVAR',
        ]
        for line in data:
            f.write(line + '\n')

@Project.pre(created_dirs)
@Project.post(built_proteins)
@Project.post(packed_proteins)
@Project.post(built_system_input)
@Project.operation
def build_tleap_input(job):
    with open(job.path + '/0_build/tleap_protein_1.in', 'w') as f:
        data = [
            'source leaprc.protein.ff14SBonlysc', 
            f'mol = sequence {{{single_to_triple(job.sp.sequences[0])}}}',
            'savepdb mol seq1.pdb',
            'quit',
        ]
        for line in data:
            f.write(line + '\n')
    
    with open(job.path + '/0_build/tleap_protein_2.in', 'w') as f:
        data = [
            'source leaprc.protein.ff14SBonlysc', 
            f'mol = sequence {{{single_to_triple(job.sp.sequences[1])}}}',
            'savepdb mol seq2.pdb',
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
    if job.sp.sequences[0] == job.sp.sequences[1]:
        with open(job.path + '/0_build/packmol.inp', 'w') as f:
            data = [
                'tolerance 40.0',
                'filetype pdb',
                'output system.pdb',
                'structure seq1_noH.pdb',
                '   number 2',
                '   inside box 0. 0. 0. 150. 150. 150.',
                'end structure'
            ]
            for line in data:
                f.write(line + '\n')
    else:
        with open(job.path + '/0_build/packmol.inp', 'w') as f:
            data = [
                'tolerance 40.0',
                'filetype pdb',
                'output system.pdb',
                'structure seq1_noH.pdb',
                '   number 1',
                '   inside box 0. 0. 0. 150. 150. 150.',
                'end structure'
                'structure seq2_noH.pdb',
                '   number 1',
                '   inside box 0. 0. 0. 150. 150. 150.',
                'end structure',
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
    total_charge = sum(aa_charge[aa1]+aa_charge[aa2] for aa1, aa2 in zip(job.sp.sequences[0], job.sp.sequences[1]))
    if total_charge == 0:
        charge_line = ''
    else:
        charge_line = f'addIons agg Cl- {total_charge}' if total_charge > 0 else f'addIons agg Na+ {abs(total_charge)}'
    with open(job.path + '/0_build/tleap_system.in', 'w') as f:
        data = [
            'source leaprc.protein.ff14SB', 
            'source leaprc.gaff', 
            'source leaprc.water.tip3p',
            'agg = loadpdb system.pdb',
            charge_line,
            'solvatebox agg TIP3PBOX 15',
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
    return f". ~/load_qz.sh; cd 0_build; tleap -f tleap_protein_1.in; pdb4amber --nohyd -i seq1.pdb -o seq1_noH.pdb; tleap -f tleap_protein_2.in; pdb4amber --nohyd -i seq2.pdb -o seq2_noH.pdb"

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
    line_num = search_pattern_in_consecutive_lines(
        '0_build/system.pdb',
        pattern1=f'  {len(job.sp.sequences[0])}   ',
        pattern2=f'  1   ',
    )
    add_line_at_line_number(
        '0_build/system.pdb',
        line_number=line_num,
        line_content='TER'
    )

    return f". ~/load_qz.sh; cd 0_build; tleap -f tleap_system.in; cp system.prmtop system.inpcrd ..; cp solvated_system.pdb ../5_meta/"

@Project.pre(created_dirs)
@Project.pre(built_proteins)
@Project.pre(packed_proteins)
@Project.pre(built_system_input)
@Project.pre(built_system)
@Project.operation(cmd=True)
def run_simulation(job):
    job.open()
    return "sbatch run_simulation.sh"

if __name__ == "__main__":
    Project().main()
