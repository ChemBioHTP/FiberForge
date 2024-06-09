import os
import shutil
import pathlib
from copy import deepcopy

import numpy as np
from flow import FlowProject
from flow.environment import DefaultSlurmEnvironment

from fiberForge.characterization import calculate_variable_over_time, estimate_elastic_modulus
from fiberForge.geometry_analysis import identify_protofibrils, identify_growth_axis
from fiberForge.utils import (
    rename_amino_acid, 
    renumber_residues, 
    remove_ligands_water_ions, 
    add_hydrogens, 
    remove_chains, 
    extract_chain, 
    calculate_n_residues, 
    discontinuous_residue_position,
    find_bounding_box
)
from fiberForge.build import calculate_cross_sectional_area

class Project(FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""
    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()

class Accre(DefaultSlurmEnvironment):
    """Subclass of DefaultPBSEnvironment for VU's ACCRE cluster."""
    # template = "accre.sh"


######################
### Calculation Labels
######################
@Project.label
def ran_preprocess_pdb(job):
    return job.isfile('0_preprocess/seq.pdb')
@Project.label
def ran_preprocess_pull(job):
    return job.isfile('4_smd/pull.tpr')
@Project.label
def finished_npt_simulation(job):
    return job.isfile('3_eq_npt/npt.gro')
@Project.label
def ran_pull(job):
    return job.isfile('4_smd/pull.xtc')
@Project.label
def finished_pull(job):
    return job.isfile('4_smd/pull.gro')
@Project.label
def ran_analysis(job):
    return 'cross_sectional_area' in job.doc

######################
### Operations
######################
@Project.post.isfile('0_preprocess/solvated_ions.gro')
@Project.operation
def preprocess_pdb(job):
    """
    Preprocess the PDB file for simulation:
        1. Remove water, ions, and ligands
        2. Unbundle fibril
        3. Identify fibril axis
        4. Estimate bounding box needed for simulation
        5. Center protein in box
        6. Solvate protein
        7. Add ions

    Steps:
    1. **Check if preprocessing is already done**:
       - If the file `0_preprocess/protofibril.gro` exists, preprocessing is skipped.
    
    2. **Set up the preprocessing directory**:
       - Create the `0_preprocess` directory if it doesn't exist.
       - Copy the `ions.mdp` file to the `0_preprocess` directory.

    3. **Check for discontinuous residue positions**:
       - Ensure there are no discontinuous residue positions in the input PDB file.
    
    4. **Preprocess the PDB file**:
       - **Rename amino acids**: Rename non-standard amino acids to standard ones.
       - **Renumber residues**: Ensure residues are numbered sequentially.
       - **Remove ligands, water, and ions**: Clean the PDB file by removing unwanted molecules.
       - **Add hydrogens**: Add missing hydrogen atoms to the protein.

    5. **Unbundle fibril**:
       - Identify protofibrils based on the fibril type and specified distance threshold.
       - Save the protofibrils information in the job document.
       - Remove chains of other protofibrils, keeping only one protofibril.

    6. **Atomtype the protein**:
       - Use `gmx pdb2gmx` to generate a GROMACS-compatible GRO file from the PDB file.

    7. **Center fibril in the box**:
       - Use `gmx editconf` to center the fibril in the simulation box.

    8. **Align fibril along the z-axis** (if required):
       - Optional step to align the fibril along the z-axis using the principal axes.

    9. **Estimate bounding box needed for simulation**:
       - Calculate the minimum and maximum coordinates to find the bounding box.
       - Identify the growth axis of the fibril.
       - Determine the pulling length and box size required for the simulation.

    10. **Expand the box to the required size**:
        - Use `gmx editconf` to expand the simulation box and center the fibril.

    11. **Solvate the protein**:
        - Use `gmx solvate` to add water molecules around the protein in the box.

    12. **Add ions**:
        - Prepare the system with `gmx grompp`.
        - Use `gmx genion` to add ions and neutralize the system.

    Exception Handling:
    - Any errors encountered during preprocessing will be caught, and a failure message will be printed.
    """
    try:
        if job.isfile(f'0_preprocess/protofibril.gro'):
            print(f"PDB file for {job.id} already preprocessed")
            return

        print(f"Preprocessing PDB file for {job.id}")

        if not os.path.exists(job.path + '/0_preprocess'):
            os.mkdir(job.path + '/0_preprocess')
        shutil.copy(job.path + "/../../in_files/ions.mdp", job.path + '/0_preprocess')
        os.chdir(job.path + '/0_preprocess')

        # Validity checks for the PDB file
        if discontinuous_residue_position(f"{job.sp.fiberverse_directory}/{job.sp.pdbID}.pdb"):
            print(f"Discontinuous residue position detected in {job.sp.pdbID}.pdb")
            return

        # Preprocess the PDB file
        input_pdb = f"{job.sp.fiberverse_directory}/{job.sp.pdbID}.pdb"
        renamed_pdb = f"{job.sp.pdbID}_renamed.pdb"
        renumbered_pdb = f"{job.sp.pdbID}_renumbered.pdb"
        cleaned_pdb = f"{job.sp.pdbID}_cleaned.pdb"
        output_pdb = f"{job.sp.pdbID}_processed.pdb"
        rename_amino_acid(input_pdb, renamed_pdb)
        renumber_residues(renamed_pdb, renumbered_pdb)
        remove_ligands_water_ions(renumbered_pdb, cleaned_pdb)
        add_hydrogens(cleaned_pdb, output_pdb)


        # Unbundle fibril
        if job.sp.type == "single":
            protofibrils = identify_protofibrils(output_pdb, distance_threshold=job.sp.protofibril_distance_threshold)
        elif job.sp.type == "solenoid":
            protofibrils = identify_protofibrils(output_pdb, distance_threshold=30.0)
        else:
            print(f"{job.sp.type} fibril type not supported yet")
            return

        # save single protofibril
        job.doc['protofibrils'] = protofibrils

        if len(protofibrils) > 1:
            chains_to_remove = []
            for i in range(1, len(protofibrils)):
                chains_to_remove.extend([chain_id for chain_id in protofibrils[i].keys()])
            remove_chains(output_pdb, "protofibril.pdb", chains_to_remove)
        else:
            os.rename(output_pdb, "protofibril.pdb")

        print("Finished identification of protofibrils")
        
        # Atomtype the protein
        os.system(f". ~/load_gromacs.sh; gmx pdb2gmx -f protofibril.pdb -o protofibril.gro -ff amber03 -water tip3p -ignh")

        # Center fibril in box
        os.system(f". ~/load_gromacs.sh; gmx editconf -f protofibril.gro -o centered.gro -c")

        align_fibril = False #TODO  # Align fibril along the z-axis
        if align_fibril:
            # Align fibril along the z-axis, assuming longest axis is the growth axis
            os.system(f". ~/load_gromacs.sh; gmx editconf -f centered.gro -o aligned.gro -princ 1 <<EOF\n1\nEOF")
        else:
            os.system(f"cp centered.gro aligned.gro")

        # Estimate bounding box needed for simulation
        min_coords, max_coords = find_bounding_box("aligned.gro")
        job.doc['min_coords'] = min_coords
        job.doc['max_coords'] = max_coords

        # Estimate fibril growth axis
        sheet_centers = list(protofibrils[0].values())
        growth_axis = identify_growth_axis(sheet_centers)
        job.doc['growth_axis'] = growth_axis

        # Calculate fibril dimensions
        fiber_container = (np.array(max_coords) - np.array(min_coords))

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
        os.system(f". ~/load_gromacs.sh; gmx editconf -f aligned.gro -o expanded.gro -box {box[0]} {box[1]} {box[2]} -center {new_center[0]} {new_center[1]} {new_center[2]}")

        # Solvate the protein
        os.system(f". ~/load_gromacs.sh; gmx solvate -cp expanded.gro -cs spc216.gro -p topol.top -o solvated.gro")

        # Add ions
        os.system(f". ~/load_gromacs.sh; gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr")
        os.system(f". ~/load_gromacs.sh; echo SOL | gmx genion -s ions.tpr -o solvated_ions.gro -p topol.top -pname NA -nname CL -neutral")

    except:
        print(f"Preprocessing failed for {job.id}")
        return

@Project.pre.isfile('0_preprocess/solvated_ions.gro')
@Project.post.isfile('run_equilibration.sh')
@Project.operation
def create_eq_submission(job):
    """
    Create and prepare a SLURM submission script for equilibration simulation.

    This function performs the following tasks:
    1. Check if the equilibration submission script already exists.
    2. Create directories for different equilibration stages.
    3. Generate the SLURM job submission script `run_equilibration.sh` which includes:
        - Job configuration details
        - Loading necessary modules and environments
        - Copying MDP files for minimization, NVT, and NPT equilibration
        - Running GROMACS commands for energy minimization, NVT, and NPT equilibration.

    Steps:
    1. **Check for existing script**:
        - If `run_equilibration.sh` already exists, the function prints a message and exits.

    2. **Create necessary directories**:
        - Creates directories `1_min`, `2_eq_nvt`, `3_eq_npt`, and `4_smd` within the job path.

    3. **Write the SLURM submission script**:
        - The script sets up the job configuration for the SLURM scheduler, including job name, partition, GPU usage, memory allocation, and time limit.
        - It loads the required modules and activates the conda environment.
        - It copies the minimization and equilibration MDP files to the respective directories.
        - It prepares and runs GROMACS commands for each equilibration stage: energy minimization, NVT equilibration, and NPT equilibration.

    Exception Handling:
    - The function does not include explicit exception handling, relying on the assumption that required directories and files are accessible and the environment is properly configured.
    """
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

@Project.pre.isfile('0_preprocess/solvated_ions.gro')
@Project.post.isfile('3_eq_npt/npt.gro')
@Project.operation(cmd=True)
def run_equilibration(job):
    job.open()
    return "sbatch run_equilibration.sh"

@Project.pre.isfile('0_preprocess/solvated_ions.gro')
@Project.pre.isfile('3_eq_npt/npt.gro')
@Project.post.isfile('4_smd/pull.mdp')
@Project.operation
def create_pull_mdp(job):
    """
    Create a GROMACS pull MDP file for umbrella sampling of a pulling simulation.

    This function performs the following tasks:
    1. Determine the reference atoms for the pulling groups.
    2. Generate the pull MDP file with appropriate parameters for umbrella sampling.

    Steps:
    1. **Determine reference atoms**:
        - Extracts chains A and B from the protofibril PDB file.
        - Calculates the center of mass (COM) for each chain.
        - Finds the closest atom to the COM for each chain, which will serve as the reference atom.

    2. **Generate pull MDP file**:
        - Constructs the pull MDP file with the following parameters:
            - Title and basic run parameters.
            - Output frequency and file settings.
            - Bond parameters and cutoff scheme.
            - PME electrostatics parameters.
            - Temperature and pressure coupling settings.
            - Pulling parameters for umbrella sampling:
                - Defines two groups (Chain A and Chain B) for pulling.
                - Specifies the reference atoms and geometry for pulling.
                - Sets the pulling rate and force constant.

    Exception Handling:
    - The function does not include explicit exception handling as it assumes necessary input files are present and calculations are valid.
    """
    # find atom that is closest to the center of mass of the pull group and use that as the reference atom
    chain_a_name = list(job.doc['protofibrils'][0].keys())[job.sp.pull_chains[0]]
    chain_b_name = list(job.doc['protofibrils'][0].keys())[job.sp.pull_chains[1]]
    chain_a = extract_chain(job.fn('0_preprocess/protofibril.pdb'), chain_a_name)
    chain_b = extract_chain(job.fn('0_preprocess/protofibril.pdb'), chain_b_name)
    chain_a_com = job.doc['protofibrils'][0][chain_a_name]
    chain_b_com = job.doc['protofibrils'][0][chain_b_name]
    # find closest atom to COM
    def find_closest_atom(com, chain):
        min_dist = np.inf
        index = 1
        count = 1
        for res in chain:
            for atom in res:
                dist = np.linalg.norm(atom.get_coord() - com)
                if dist < min_dist:
                    min_dist = dist
                    index = count
                count += 1
        return index
    ref_atom_a_index = find_closest_atom(chain_a_com, chain_a)
    ref_atom_b_index = find_closest_atom(chain_b_com, chain_b)
    n_atom_per_chain = sum(1 for atom in chain_a.get_atoms())
    assert sum(1 for atom in chain_a.get_atoms()) == sum(1 for atom in chain_b.get_atoms()) # make sure they are the same number of atoms

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
    pull_group1_name        = Chain_A 
    pull_group2_name        = Chain_B 
    pull-pbc-ref-prev-step-com = yes    ; use the reference group from the previous step since we are pulling a big group
    pull-group1-pbcatom     = {ref_atom_a_index + n_atom_per_chain * (job.sp.pull_chains[0])}
    pull-group2-pbcatom     = {ref_atom_b_index + n_atom_per_chain * (list(range(len(job.doc['protofibrils'][0])))[job.sp.pull_chains[1]])}
    pull_coord1_type        = umbrella  ; harmonic potential
    pull_coord1_geometry    = distance  ; simple distance increase 
    pull_coord1_dim         = Y Y Y  ; calculate the distance between the centers of mass in all three dimensions
    pull_coord1_groups      = 1 2
    pull_coord1_start       = yes       ; define initial COM distance > 0
    pull_coord1_rate        = {job.sp.pull_rate}      ;  nm per ps
    pull_coord1_k           = {job.sp.pull_constant}      ; kJ mol^-1 nm^-2
    """
    with open(job.path + '/4_smd/pull.mdp', 'w') as f:
        f.write(lines)

@Project.pre.isfile('3_eq_npt/npt.gro')
@Project.post.isfile('run_pullings.sh')
@Project.operation
def create_pull_submission(job):
    """
    Create a SLURM submission script for running pulling simulations.

    This function performs the following tasks:
    1. Generate a SLURM submission script (`run_pulling.sh`) for running the pulling simulation.
    2. Configure the SLURM job with appropriate parameters.

    Steps:
    1. **Generate SLURM script**:
        - Creates a SLURM submission script (`run_pulling.sh`) in the job directory.
        - The script sets up the necessary environment and modules, including AMBER, GCC, OpenMPI, and GROMACS.
        - Activates the `fiberForge` conda environment.
        - Directs to the directory containing the pull MDP file (`4_smd`).
        - Executes the GROMACS `mdrun` command to perform the pulling simulation.
        - Prints a message indicating the completion of the pulling simulation.

    Exception Handling:
    - The function assumes the availability of required input files and properly configured environment.
    """
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

@Project.pre.isfile('3_eq_npt/npt.gro')
@Project.pre.isfile('4_smd/pull.mdp')
@Project.post.isfile('4_smd/pull.tpr')
@Project.operation
def preprocess_pull(job):
    """
    Preprocesses the system for the pulling simulation.

    This function performs the following tasks:
    1. Removes positional restraints for chains except the chain to be pulled.
    2. Generates an index file for the pulling simulation.
    3. Creates the GROMACS input files for the pulling simulation.

    Steps:
    1. **Check prerequisite files**:
        - Verifies the availability of the equilibrium NPT structure (`3_eq_npt/npt.gro`) and the pull MDP file (`4_smd/pull.mdp`).

    2. **Remove positional restraints**:
        - Removes positional restraints for all chains except the chain intended for pulling. This step prepares the system for pulling.

    3. **Generate index file**:
        - Executes `make_ndx` GROMACS command to generate an index file (`index.ndx`) based on the equilibrium structure.
        - Defines two groups corresponding to the chains involved in the pulling simulation.

    4. **Create pull code**:
        - Uses `grompp` to create the GROMACS portable binary run file (`pull.tpr`) for the pulling simulation.
        - Incorporates the equilibrium structure, topology, index file, and pull MDP settings into the binary file.

    Exception Handling:
    - If prerequisite files are missing, the function prints a message indicating that the job is not ready for pull preprocessing and returns without further processing.
    - No explicit error handling is implemented beyond this point, assuming successful completion of individual steps.
    """
    if not job.isfile(f'3_eq_npt/npt.gro'):
        print("Job not ready for pull preprocessing")
        return
    os.chdir(job.path + '/4_smd')

    # Calculate the number of residues in a chain of the protofibril
    n_chains = len(job.doc['protofibrils'][0]) # Number of chains in the first protofibril
    chain = extract_chain('../0_preprocess/protofibril.pdb', list(job.doc['protofibrils'][0].keys())[0])
    n_residues = calculate_n_residues(chain)
    chain_b = list(range(n_chains))[job.sp.pull_chains[1]] + 1 # 1-indexed chain number 
    
    # Remove positional restraints for all chains restraint chain
    for chain in job.doc['protofibrils'][0].keys():
        if chain != list(job.doc['protofibrils'][0].keys())[job.sp.pull_chains[1]]:
            os.system(f"sed -i '/; Include Position restraint file/,/#endif/d' ../0_preprocess/topol_Protein_chain_{chain}.itp")

    # Make index file
    ## Run make_ndx once to figure out how many groups there are
    os.system(f". ~/load_gromacs.sh; gmx make_ndx -f ../3_eq_npt/npt.gro -o index.ndx<<EOF\nq\nEOF")
    ## Process index.ndx to find the number of groups
    n_groups = 0
    with open('index.ndx', 'r') as f:
        for line in f:
            if "[" in line:
                n_groups += 1
    os.system(f". ~/load_gromacs.sh; gmx make_ndx -f ../3_eq_npt/npt.gro -o index.ndx<<EOF\nri 1-{n_residues}\nname {n_groups} Chain_A\nri {n_residues*(chain_b - 1)+1}-{n_residues * chain_b}\nname {n_groups+1} Chain_B\nq\nEOF")

    # Create the pull code
    os.system(". ~/load_gromacs.sh; gmx grompp -f pull.mdp -c ../3_eq_npt/npt.gro -p ../0_preprocess/topol.top -r ../3_eq_npt/npt.gro -n index.ndx -t ../3_eq_npt/npt.cpt -o pull.tpr")

@Project.pre.isfile('run_pulling.sh')
@Project.post(ran_pull)
@Project.operation(cmd=True)
def run_pull(job):
    job.open()
    return "sbatch run_pulling.sh"

@Project.pre(ran_pull)
@Project.post(lambda job: job.doc.get('stress'))
@Project.operation
def run_analysis(job):
    """
    Runs analysis on the results of the pulling simulation.

    This function performs the following analyses:
    1. Calculates the cross-sectional area of the fibril.
    2. Determines the strain and stress over time during pulling.
    3. Computes the ultimate tensile strength and elastic modulus.

    Steps:
    1. **Calculate Cross-sectional Area**:
        - Computes the cross-sectional area of the fibril perpendicular to the growth axis.
        - Uses the `calculate_cross_sectional_area` function with the fibril axis and the protofibril PDB file.

    2. **Determine Strain and Stress**:
        - Calculates the strain over time from the pulling trajectory (`pull_pullx.xvg`).
        - Converts the force over time (`pull_pullf.xvg`) to stress using the cross-sectional area.
        - Stores strain and stress data in the job document.

    3. **Compute Ultimate Tensile Strength and Elastic Modulus**:
        - Finds the maximum stress to determine the ultimate tensile strength.
        - Calculates the elastic modulus by dividing the maximum stress by the strain at that point.
        - Stores the computed values in the job document.

    Exception Handling:
    - If any step fails, the function prints a message indicating the failure and returns without further processing.
    """
    # Calculate the cross-sectional area
    if not os.path.exists(job.path + '/0_preprocess/protofibril.pdb'):
        print(f"Protofibril.pdb not found for {job.id}")
        return
    else:
        print(f"Running analysis on {job.id}")
        fibril_axis = [0, 0, 0]
        fibril_axis[job.doc['growth_axis']] = 1
        cross_sectional_area = calculate_cross_sectional_area(fibril_axis, job.path+ '/0_preprocess/protofibril.pdb')
        cross_sectional_area = cross_sectional_area * 1e-20 # A^2 to m^2
        job.doc['cross_sectional_area'] = cross_sectional_area

        # Calculate the strain in units of nm
        time_length = calculate_variable_over_time(job.path + '/4_smd/pull_pullx.xvg')
        length_over_time = np.array([l for (t, l) in time_length])
        strain = (length_over_time - length_over_time[0]) / length_over_time[0]
        job.doc['strain'] = strain # nm/nm of pull distance

        # Calculate the stress
        time_force = calculate_variable_over_time(job.path + '/4_smd/pull_pullf.xvg')
        force_over_time = np.array([f for (t, f) in time_force]) * (1e-9) * (1/1000) # kJ/mol/nm to N
        stress = force_over_time / cross_sectional_area
        job.doc['stress'] = stress

        # Calculate the ultimate tensile strength
        ultimate_tensile_strength = np.max(stress)
        job.doc['ultimate_tensile_strength'] = ultimate_tensile_strength

        # Calculate the elastic modulus, assuming there is no plastic deformation
        E, yield_point = estimate_elastic_modulus(stress, strain)
        job.doc['elastic_modulus'] = E
        job.doc['yield_point'] = yield_point

@Project.pre.isfile('4_smd/pull.trr')
@Project.post(lambda job: not job.isfile('4_smd/pull.trr'))
@Project.operation
def remove_trr(job):
    try:
        os.remove(job.fn('4_smd/pull.trr'))
        os.remove(job.fn('3_eq_npt/npt.trr'))
        os.remove(job.fn('2_eq_nvt/nvt.trr'))
    except:
        pass

if __name__ == "__main__":
    Project().main()
