import os
import shutil
import pathlib
from copy import deepcopy

import numpy as np
from flow import FlowProject
from flow.environment import DefaultSlurmEnvironment
from textwrap import wrap
import signac
from sklearn.decomposition import PCA
import pymol2
from Bio.PDB import PDBParser, PDBIO, Select
from fiberForge.utils import *

import warnings
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



def remove_ligands_water_ions(input_pdb, output_file):
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(input_pdb)

        # Remove ligands
        try:
            pymol.cmd.remove('het')
        except:
            print('No ligands found')
            
        # Remove water
        try:
            pymol.cmd.remove('solvent')
        except:
            print('No water found')

        # Remove ions
        try:
            pymol.cmd.remove('ions')
        except:
            print('No ions found')

        # Save the modified structure
        pymol.cmd.save(output_file)

def add_hydrogens(input_pdb, output_pdb):
    with pymol2.PyMOL() as pymol:
        # Load the input PDB file
        pymol.cmd.load(input_pdb)

        # Add hydrogens
        pymol.cmd.h_add()

        # Save the modified structure
        pymol.cmd.save(output_pdb)

def find_bounding_box(gro_file):
    with open(gro_file, 'r') as f:
        lines = f.readlines()
        n_atoms = int(lines[1])
        coords = np.zeros((n_atoms, 3))
        for i in range(2, n_atoms + 1):
            coords[i - 2] = np.array([float(x) for x in lines[i].split()[-3:]])
    min_coords = np.min(coords, axis=0)
    max_coords = np.max(coords, axis=0)
    return min_coords, max_coords

def identify_protofibrils(pdb_file, distance_threshold=20.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    calculate_distance = lambda c1, c2: np.linalg.norm(c1 - c2)
    calculate_center_of_mass = lambda coords: np.mean(coords, axis=0)

    chain_centers = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id
            coords = []
            for residue in chain:
                for atom in residue:
                    coords.append(atom.get_coord())
            coords = np.array(coords)
            center_of_mass = calculate_center_of_mass(coords)
            chain_centers[chain_id] = center_of_mass

    protofibrils = []

    for chain_id, center in chain_centers.items():
        found = False
        for protofibril in protofibrils:
            for chain_id_in_protofibril, center_in_protofibril in protofibril.items():
                distance = calculate_distance(center, center_in_protofibril)
                if distance < distance_threshold:
                    protofibril[chain_id] = center
                    found = True
                    break
            if found:
                break

        if not found:
            protofibrils.append({chain_id: center})

    return protofibrils

def remove_chains(input_pdb, output_pdb, chains_to_keep):
    class ChainSelector(Select):
        def accept_chain(self, chain):
            return chain.id in chains_to_keep

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=ChainSelector())

def calculate_fibril_dimensions(coords):
    pca = PCA(n_components=3)
    pca.fit(coords)
    return pca.explained_variance_

def extract_chain(input_pdb, chain_id):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)

    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                return chain

    return None

def identify_growth_axis(sheet_centers):
    # Convert sheet centers to a NumPy array
    sheet_centers = np.array(sheet_centers)

    # Calculate the covariance matrix
    covariance_matrix = np.cov(sheet_centers, rowvar=False)

    # Calculate the eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)

    # Identify the index of the maximum eigenvalue (corresponding to the growth axis)
    growth_axis = np.argmax(eigenvalues)

    return growth_axis

def calculate_n_residues(chain):
    n_residues = 0
    for residue in chain:
        n_residues += 1
    return n_residues

def renumber_residues(input_pdb, output_pdb):
    class RenumberSelector(Select):
        def __init__(self):
            self.chain_id_mapping = {}  # Mapping of old chain IDs to new chain IDs
            self.residue_number_mapping = {}  # Mapping of old residue numbers to new residue numbers
            self.new_chain_id = 'A'     # Starting new chain ID
            self.new_residue_number = 1  # Starting new residue number

        def accept_chain(self, chain):
            # Update the chain ID and store the mapping
            new_chain_id = self.chain_id_mapping.get(chain.id, None)
            if new_chain_id is None:
                new_chain_id = self.new_chain_id
                self.chain_id_mapping[chain.id] = new_chain_id
                self.new_chain_id = chr(ord(self.new_chain_id) + 1)

            # Update the residue numbers and store the mapping
            for residue in chain:
                new_residue_number = self.residue_number_mapping.get(residue.id[1], None)
                if new_residue_number is None:
                    new_residue_number = self.new_residue_number
                    self.residue_number_mapping[residue.id[1]] = new_residue_number
                    self.new_residue_number += 1

                residue.id = (' ', new_residue_number, ' ')

            return True

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=RenumberSelector())

def discontinuous_residue_position(pdb_file, distance_threshold=10.0):
    remove_ligands_water_ions(pdb_file, 'temp.pdb')
    structure = PDBParser(QUIET=True).get_structure('protein', 'temp.pdb')
    os.remove('temp.pdb')

    for model in structure:
        for chain in model:
            previous_position = None

            for residue in chain:
                current_position = np.array(residue['CA'].coord)

                if previous_position is not None:
                    distance = np.linalg.norm(current_position - previous_position)

                    if distance > distance_threshold:
                        print(f"Chain {chain.id}, discontinuous change detected between residues {residue.id[1]} and {previous_residue.id[1]}")
                        return True

                previous_position = current_position
                previous_residue = residue
    return False

def rename_amino_acid(pdb_file_path, output_pdb_path):
    structure = PDBParser(QUIET=True).get_structure('protein', pdb_file_path)

    for model in structure:
        for chain in model:
            for residue in chain:
                original_residue_name = residue.resname.strip()

                # Check for D-proline and rename to normal Proline
                if original_residue_name == 'DPR':
                    residue.resname = 'PRO'

    # Save the modified structure to a new PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb_path)

    print(f"Modified PDB file saved to: {output_pdb_path}")

def extract_region(input_pdb, output_pdb, region):
    class RegionSelector(Select):
        def accept_residue(self, residue):
            return region[0] <= residue.id[1] <= region[1]

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=RegionSelector())

def n_residues(pdb_file):
    structure = PDBParser(QUIET=True).get_structure('protein', pdb_file)
    n_residues = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                n_residues += 1
    return n_residues

"""
Preprocess the PDB file
    1. Remove water, ions, and ligands
    2. Unbundle fibril
    3. Identify fibril axis
    4. Estimate bounding box needed for simulation
    5. Center protein in box
    6. Solvate protein
    7. Add ions
"""

######################
### Calculation Labels
######################
# @Project.label
# def ran_preprocess_pdb(job):
#     return job.isfile('0_preprocess/seq.pdb')
# @Project.label
# def ran_preprocess_pull(job):
#     return job.isfile('4_smd/pull.tpr')
@Project.label
def finished_npt_simulation(job):
    return job.isfile('3_eq_npt/npt.gro')
@Project.label
def ran_pull(job):
    return job.isfile('4_smd/pull.xtc')
@Project.label
def finished_pull(job):
    return job.isfile('4_smd/pull.gro')
######################
### Operations
######################
@Project.post.isfile('0_preprocess/solvated_ions.gro')
@Project.operation
def preprocess_pdb(job):
    if job.isfile(f'0_preprocess/protofibril.gro'):
        print(f"PDB file for {job.id} already preprocessed")
        return
    
    print(f"Preprocessing PDB file for {job.id}")

    if not os.path.exists(job.path + '/0_preprocess'):
        os.mkdir(job.path + '/0_preprocess')
    shutil.copy(job.path + "/../../in_files/ions.mdp", job.path + '/0_preprocess')
    os.chdir(job.path + '/0_preprocess')

    # Cut the region of interest from the PDB file
    input_pdb = f"{job.sp.data_dir}/{job.sp.pdbID}.pdb"
    cut_pdb = f"{job.sp.pdbID}_cut.pdb"
    chains_pdb = f"{job.sp.pdbID}_chains.pdb"
    remove_chains(input_pdb, chains_pdb, [job.sp.chain])
    extract_region(chains_pdb, cut_pdb, job.sp.region)

    # Validity checks for the PDB file
    if discontinuous_residue_position(f"{job.sp.data_dir}/{job.sp.pdbID}.pdb"):
        print(f"Discontinuous residue position detected in {job.sp.pdbID}.pdb")
        return

    # Preprocess the PDB file
    renamed_pdb = f"{job.sp.pdbID}_renamed.pdb"
    renumbered_pdb = f"{job.sp.pdbID}_renumbered.pdb"
    cleaned_pdb = f"{job.sp.pdbID}_cleaned.pdb"
    output_pdb = f"{job.sp.pdbID}_processed.pdb"
    rename_amino_acid(cut_pdb, renamed_pdb)
    renumber_residues(renamed_pdb, renumbered_pdb)
    remove_ligands_water_ions(renumbered_pdb, cleaned_pdb)
    add_hydrogens(cleaned_pdb, output_pdb)

    # Atomtype the protein
    os.system(f". ~/load_gromacs.sh; gmx pdb2gmx -f {output_pdb} -o out.gro -ff amber03 -water tip3p -ignh")

    # Center fibril in box
    os.system(f". ~/load_gromacs.sh; gmx editconf -f out.gro -o centered.gro -c")

    # Estimate bounding box needed for simulation
    min_coords, max_coords = find_bounding_box("centered.gro")
    job.doc['min_coords'] = min_coords
    job.doc['max_coords'] = max_coords    

    # Calculate fibril dimensions
    fiber_container = (np.array(max_coords) - np.array(min_coords))

    # Estimate fibril growth axis
    growth_axis = np.argmax(fiber_container)
    job.doc['growth_axis'] = growth_axis

    # Create the pulling length, for now 5 times the length of the protofibril
    pulling_length = 5 * fiber_container[growth_axis]
    job.doc['pulling_length'] = pulling_length

    # Define the box size as the pulling length in the growth axis and 2 times the other dimensions in the other axes
    box = deepcopy(fiber_container)
    box[growth_axis] = pulling_length
    box[box != pulling_length] = 2 * box[box != pulling_length]
    job.doc['box_size'] = box

    # Find new center for fibril
    new_center = deepcopy(box)
    new_center[growth_axis] = new_center[growth_axis] / 5
    new_center[box != pulling_length] = box[box != pulling_length] / 2

    # Expand the box to the box size
    os.system(f". ~/load_gromacs.sh; gmx editconf -f centered.gro -o expanded.gro -box {box[0]} {box[1]} {box[2]} -center {new_center[0]} {new_center[1]} {new_center[2]}")

    # Solvate the protein
    os.system(f". ~/load_gromacs.sh; gmx solvate -cp expanded.gro -cs spc216.gro -p topol.top -o solvated.gro")

    # Add ions
    os.system(f". ~/load_gromacs.sh; gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr")
    os.system(f". ~/load_gromacs.sh; echo SOL | gmx genion -s ions.tpr -o solvated_ions.gro -p topol.top -pname NA -nname CL -neutral")

    # Remove positional restraints for all chains except the second chain

@Project.pre.isfile('3_eq_npt/npt.gro')
@Project.pre.isfile('4_smd/pull.mdp')
@Project.post.isfile('4_smd/pull.tpr')
@Project.operation
def preprocess_pull(job):
    if not job.isfile(f'3_eq_npt/npt.gro'):
        print("Job not ready for pull preprocessing")
        return
    os.chdir(job.path + '/4_smd')

    # Make index file
    # Run make_ndx once to figure out how many groups there are
    os.system(f". ~/load_gromacs.sh; gmx make_ndx -f ../3_eq_npt/npt.gro -o index.ndx<<EOF\nq\nEOF")
    
    # Process index.ndx to find the number of groups
    n_groups = 0
    with open('index.ndx', 'r') as f:
        for line in f:
            if "[" in line:
                n_groups += 1
    n_residues_pg_1 = job.sp.pull_group1[1] - job.sp.pull_group1[0]
    n_residues_pg_2 = job.sp.pull_group2[1] - job.sp.pull_group2[0]
    pg_2_start = n_residues(f"../0_preprocess/{job.sp.pdbID}_processed.pdb") + job.sp.pull_group2[0]
    pg_2_end = pg_2_start + n_residues_pg_2
    os.system(f". ~/load_gromacs.sh; gmx make_ndx -f ../3_eq_npt/npt.gro -o index.ndx<<EOF\nri 1-{n_residues_pg_1}\nname {n_groups} Group_A\nri {pg_2_start}-{pg_2_end}\nname {n_groups+1} Group_B\nq\nEOF")

    # Create the position restraint file
    os.system(f"gmx genrestr -f ../3_eq_npt/npt.gro -n index.ndx -o ../0_preprocess/posre.itp -fc 1000 1000 1000 <<EOF\n{n_groups}\nEOF")

    # Create the pull code
    os.system(". ~/load_gromacs.sh; gmx grompp -f pull.mdp -c ../3_eq_npt/npt.gro -p ../0_preprocess/topol.top -r ../3_eq_npt/npt.gro -n index.ndx -t ../3_eq_npt/npt.cpt -o pull.tpr")

# @Project.operation
# def create_eq_mdp(job):
#     pass

@Project.pre.isfile('0_preprocess/solvated_ions.gro')
@Project.pre.isfile('3_eq_npt/npt.gro')
@Project.post.isfile('4_smd/pull.mdp')
@Project.operation
def create_pull_mdp(job):
    lines = f"""title       = Umbrella pulling simulation 
    define      = -DPOSRES
    ; Run parameters
    integrator  = md
    dt          = 0.002
    tinit       = 0
    nsteps      = {job.sp.pull_steps}    ; 2 femtoseconds per step
    nstcomm     = 10
    ; Output parameters
    nstxout     = 500      ; every 1 ps
    nstvout     = 5000 
    nstfout     = 500
    nstxtcout   = 500       ; every 1 ps
    nstenergy   = 500
    ; Bond parameters
    constraint_algorithm    = lincs
    constraints             = all-bonds
    continuation            = yes       ; continuing from NPT 
    ; Single-range cutoff scheme
    cutoff-scheme   = Verlet
    nstlist         = 20 
    ns_type         = grid 
    rlist           = 1.4
    rcoulomb        = 1.4
    rvdw            = 1.4
    ; PME electrostatics parameters
    coulombtype     = PME
    fourierspacing  = 0.12
    fourier_nx      = 0
    fourier_ny      = 0
    fourier_nz      = 0
    pme_order       = 4
    ewald_rtol      = 1e-5
    optimize_fft    = yes
    ; Berendsen temperature coupling is on in two groups
    Tcoupl      = Nose-Hoover
    tc_grps     = Protein   Non-Protein 
    tau_t       = 1.0       1.0
    ref_t       = 310       310
    ; Pressure coupling is on
    Pcoupl          = Parrinello-Rahman 
    pcoupltype      = isotropic
    tau_p           = 2.0       
    compressibility = 4.5e-5
    ref_p           = 1.0
    refcoord_scaling = com
    ; Generate velocities is off
    gen_vel     = no 
    ; Periodic boundary conditions are on in all directions
    pbc     = xyz
    ; Long-range dispersion correction
    DispCorr    = EnerPres
    ; Pull code
    pull                    = yes
    pull_ncoords            = 1         ; only one reaction coordinate 
    pull_ngroups            = 2         ; two groups defining one reaction coordinate 
    pull_group1_name        = Group_A 
    pull_group2_name        = Group_B 
    pull_coord1_type        = umbrella  ; harmonic potential
    pull_coord1_geometry    = distance  ; simple distance increase 
    pull_coord1_dim         = N N Y
    pull_coord1_groups      = 1 2
    pull_coord1_start       = yes       ; define initial COM distance > 0
    pull_coord1_rate        = {job.sp.pull_rate}      ;  nm per ps
    pull_coord1_k           = {job.sp.pull_constant}      ; kJ mol^-1 nm^-2
    """
    with open(job.path + '/4_smd/pull.mdp', 'w') as f:
        f.write(lines)

@Project.pre.isfile('0_preprocess/solvated_ions.gro')
@Project.post.isfile('run_equilibration.sh')
@Project.operation
def create_eq_submission(job):
    if os.path.exists(job.path + '/run_equilibration.sh'):
        print("Equilibration submission script already exists")
        return
    os.mkdir(job.path + '/1_min')
    os.mkdir(job.path + '/2_eq_nvt')
    os.mkdir(job.path + '/3_eq_npt')
    os.mkdir(job.path + '/4_smd')
    with open(job.path + '/run_equilibration.sh', 'w') as f:
        data = [
            "#!/bin/bash",
            f"#SBATCH --job-name=eq_{job.id}",
            "#SBATCH --account=csb_gpu_acc",
            "#SBATCH --partition=pascal",  
            "#SBATCH --gres=gpu:1",
            "#SBATCH --no-requeue",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks-per-node=1",
            "#SBATCH --time=72:00:00",
            "#SBATCH --mem=10G",            
            "#SBATCH --output=eq.log",
            ""
            "source /home/shaoq1/bin/amber_env/amber-accre.sh",
            "module load GCC/8.2.0",
            "module load OpenMPI/3.1.4",
            "module load GROMACS/2020",
            "source ~/.bashrc",
            "conda activate fiberForge",
            "",
            "cp ../../in_files/minim.mdp 1_min",
            "cp ../../in_files/nvt.mdp 2_eq_nvt",
            "cp ../../in_files/npt.mdp 3_eq_npt",
            "",
            "cd 1_min",
            "gmx grompp -f minim.mdp -c ../0_preprocess/solvated_ions.gro -p ../0_preprocess/topol.top -o em.tpr",
            "gmx mdrun -v -deffnm em",
            "cd ..",
            "cd 2_eq_nvt",
            "gmx grompp -f nvt.mdp -c ../1_min/em.gro -p ../0_preprocess/topol.top -r ../1_min/em.gro -o nvt.tpr",
            "gmx mdrun -deffnm nvt",
            "cd ..",
            "cd 3_eq_npt",
            "gmx grompp -f npt.mdp -c ../2_eq_nvt/nvt.gro -p ../0_preprocess/topol.top -r ../2_eq_nvt/nvt.gro -o npt.tpr",
            "gmx mdrun -deffnm npt",
            "cd ..",
            "echo 'Finished equilibration'",
        ]
        for line in data:
            f.write(line + '\n')

@Project.pre.isfile('3_eq_npt/npt.gro')
@Project.post.isfile('run_pullings.sh')
@Project.operation
def create_pull_submission(job):
    with open(job.path + '/run_pulling.sh', 'w') as f:
        data = [
            "#!/bin/bash",
            f"#SBATCH --job-name=smd_{job.id}",
            "#SBATCH --account=csb_gpu_acc",
            "#SBATCH --partition=pascal",  
            "#SBATCH --gres=gpu:1",
            "#SBATCH --no-requeue",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks-per-node=1",
            "#SBATCH --time=72:00:00",
            "#SBATCH --mem=10G",            
            "#SBATCH --output=pull.log",
            ""
            "source /home/shaoq1/bin/amber_env/amber-accre.sh",
            "module load GCC/8.2.0",
            "module load OpenMPI/3.1.4",
            "module load GROMACS/2020",
            "source ~/.bashrc",
            "conda activate fiberForge",
            "",
            "cd 4_smd",
            "gmx mdrun -deffnm pull",
            "echo 'Finished pulling'",
        ]
        for line in data:
            f.write(line + '\n')

@Project.pre.isfile('0_preprocess/solvated_ions.gro'    )
@Project.post.isfile('3_eq_npt/npt.gro')
@Project.operation(cmd=True)
def run_equilibration(job):
    job.open()
    return "sbatch run_equilibration.sh"

@Project.pre.isfile('run_pulling.sh')
@Project.post(ran_pull)
@Project.operation(cmd=True)
def run_pull(job):
    job.open()
    return "sbatch run_pulling.sh"

if __name__ == "__main__":
    Project().main()
