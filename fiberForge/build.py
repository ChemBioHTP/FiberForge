import numpy as np
import shutil
from plotly import graph_objects as go
from scipy.optimize import minimize
import mdtraj
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue, Atom
from scipy.optimize import minimize
from sklearn.decomposition import PCA
import os
import subprocess

# def visualize_rotation_translation(pdb_file, chain1_id, chain2_id):
#     rotation, translation = estimate_rotation_translation_between_chains(pdb_file, chain1_id, chain2_id)
    
#     parser = PDBParser()
#     structure = parser.get_structure("protein", pdb_file)
    
#     fig = go.Figure()
    
#     # Extract coordinates of atoms from the two specified chains
#     coords_chain1 = []
#     coords_chain2 = []
#     for model in structure:
#         for chain in model:
#             if chain.id == chain1_id:
#                 for residue in chain:
#                     for atom in residue:
#                         coords_chain1.append(atom.get_coord())
#             elif chain.id == chain2_id:
#                 for residue in chain:
#                     for atom in residue:
#                         coords_chain2.append(atom.get_coord())
    
#     coords_chain1 = np.array(coords_chain1)
#     coords_chain2 = np.array(coords_chain2)
    
#     # Apply rotation and translation to chain1
#     transformed_coords_chain1 = np.dot(coords_chain1, rotation.T) + translation
    
#     # Plot the original and transformed coordinates
#     x_chain2, y_chain2, z_chain2 = coords_chain2[:,0], coords_chain2[:,1], coords_chain2[:,2]
#     x_transformed_chain1, y_transformed_chain1, z_transformed_chain1 = transformed_coords_chain1[:,0], transformed_coords_chain1[:,1], transformed_coords_chain1[:,2]
    
#     # fig.add_trace(go.Scatter3d(x=x_chain1, y=y_chain1, z=z_chain1, mode='markers', name='Chain 1 (Original)'))
#     fig.add_trace(go.Scatter3d(x=x_chain2, y=y_chain2, z=z_chain2, mode='markers', name='Chain 2'))
#     fig.add_trace(go.Scatter3d(x=x_transformed_chain1, y=y_transformed_chain1, z=z_transformed_chain1, mode='markers', name='Chain 1 (Transformed)'))
    
#     fig.update_layout(title='Rotation and Translation Visualization',
#                       scene=dict(
#                           xaxis_title='X',
#                           yaxis_title='Y',
#                           zaxis_title='Z'
#                       ))
    
#     fig.show()

def calculate_rmsd(chain1_coords, chain2_coords):
    """
    Calculate the Root Mean Square Deviation (RMSD) between two sets of coordinates.

    Parameters:
    - chain1_coords (np.array): Coordinates of atoms in the original chain.
    - chain2_coords (np.array): Coordinates of atoms in the second chain.

    Returns:
    - rmsd (float): The RMSD value.
    """
    # Calculate squared distances between corresponding atoms
    squared_distances = np.sum((chain1_coords - chain2_coords) ** 2, axis=1)
    # Calculate RMSD
    rmsd = np.sqrt(np.mean(squared_distances))
    return rmsd

def calculate_rmsd_between_chains(pdb_file, chain1_id, chain2_id, rotation, translation):
    """
    Calculate the RMSD between the second chain and the transformed first chain.

    Parameters:
    - pdb_file (str): Path to the PDB file.
    - chain1_id (str): ID of the first chain.
    - chain2_id (str): ID of the second chain.
    - rotation (np.array): Rotation matrix.
    - translation (np.array): Translation vector.

    Returns:
    - rmsd_transformed (float): RMSD between the transformed first chain and the second chain.
    """
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    
    # Extract coordinates of atoms from the two specified chains
    coords_chain1 = []
    coords_chain2 = []
    for model in structure:
        for chain in model:
            if chain.id == chain1_id:
                for residue in chain:
                    for atom in residue:
                        coords_chain1.append(atom.get_coord())
            elif chain.id == chain2_id:
                for residue in chain:
                    for atom in residue:
                        coords_chain2.append(atom.get_coord())
    
    coords_chain1 = np.array(coords_chain1)
    coords_chain2 = np.array(coords_chain2)
    
    # Apply rotation and translation to chain1
    transformed_coords_chain1 = np.dot(coords_chain1, rotation.T) + translation
    
    # Calculate RMSD between the transformed chain1 and chain2
    rmsd_transformed = calculate_rmsd(transformed_coords_chain1, coords_chain2)
    
    return rmsd_transformed

def calculate_residuals(chain1_coords, chain2_coords):
    """
    Calculate the residuals between two sets of coordinates.

    Parameters:
    - chain1_coords (np.array): Coordinates of atoms in the first chain.
    - chain2_coords (np.array): Coordinates of atoms in the second chain.

    Returns:
    - residuals (np.array): Array of residuals.
    """
    residuals = chain1_coords - chain2_coords
    return residuals

def visualize_residuals(pdb_file, chain1_id, chain2_id, rotation, translation):
    """
    Visualize the residuals between the second chain and the transformed first chain using 3D cones.

    Parameters:
    - pdb_file (str): Path to the PDB file.
    - chain1_id (str): ID of the first chain.
    - chain2_id (str): ID of the second chain.
    - rotation (np.array): Rotation matrix.
    - translation (np.array): Translation vector.
    """
    import random
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    
    # Extract coordinates of atoms from the two specified chains
    coords_chain1 = []
    coords_chain2 = []
    for model in structure:
        for chain in model:
            if chain.id == chain1_id:
                for residue in chain:
                    for atom in residue:
                        coords_chain1.append(atom.get_coord())
            elif chain.id == chain2_id:
                for residue in chain:
                    for atom in residue:
                        coords_chain2.append(atom.get_coord())
    
    coords_chain1 = np.array(coords_chain1)
    coords_chain2 = np.array(coords_chain2)
    
    # Apply rotation and translation to chain1
    transformed_coords_chain1 = np.dot(coords_chain1, rotation.T) + translation
    
    # Calculate residuals between chain2 and transformed chain1
    residuals = calculate_residuals(transformed_coords_chain1, coords_chain2)
    
    # Create figure
    fig = go.Figure()
    
    # Add 3D cones for each residual
    sample_size = 2000
    sample = random.sample(list(range(len(residuals))), sample_size)
    all_res = np.array(residuals)[sample]
    all_coords = np.array(coords_chain2)[sample]
    # fig.add_trace(go.Cone(x=[x], y=[y], z=[z], u=[u], v=[v], w=[w], colorscale='Blues', sizemode="absolute"))

    fig.add_trace(go.Cone(
        x=all_coords[:,0],
        y=all_coords[:,1],
        z=all_coords[:,2],
        u=all_res[:,0],
        v=all_res[:,1],
        w=all_res[:,2],
        colorscale='Jet',
        sizemode="absolute",
    ))

    fig.update_layout(
        autosize=False,
        width=500,
        height=500,
        margin=dict(
            l=50,
            r=50,
            b=100,
            t=100,
            pad=4
        ),
    )
    
    # # Add scatter plot of chain2 atoms
    # fig.add_trace(go.Scatter3d(x=coords_chain2[:,0], y=coords_chain2[:,1], z=coords_chain2[:,2], mode='markers', marker=dict(color='red'), name='Chain 2'))
    
    # Update layout
    fig.update_layout(title='Residuals between Second Chain and Transformed First Chain',
                      scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'))
    
    # Show figure
    fig.show()

def calculate_average_rotation_translation(pdb_file):
    """
    Calculate the average rotation and translation between subsequent chains in a PDB file.

    Parameters:
    - pdb_file (str): Path to the PDB file.

    Returns:
    - average_rotation (np.array): Average rotation matrix.
    - average_translation (np.array): Average translation vector.
    """
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)

    chains_coords = []
    
    # Iterate over each subsequent chain pair
    for model in structure:
        chains = list(model)
        for i in range(len(chains) - 1):
            chain1_id = chains[i].id
            chain2_id = chains[i + 1].id     
            # Extract coordinates of atoms from the two specified chains
            coords_chain1 = []
            coords_chain2 = []
            for model in structure: #TODO this is assuming the order of the chainID is the same as the order of the chains in the structure
                for chain in model:
                    if chain.id == chain1_id:
                        for residue in chain:
                            for atom in residue:
                                coords_chain1.append(atom.get_coord())
                    elif chain.id == chain2_id:
                        for residue in chain:
                            for atom in residue:
                                coords_chain2.append(atom.get_coord())
            
            coords_chain1 = np.array(coords_chain1)
            coords_chain2 = np.array(coords_chain2)

            chains_coords.append([coords_chain1, coords_chain2])

    # Define the function to minimize (sum of squared distances)
    def objective(params):
        rotation_matrix = np.reshape(params[:9], (3, 3))
        translation = params[9:]
        total_rmsd = 0.0
        for (coords_chain1, coords_chain2) in chains_coords:
            transformed_coords_chain1 = np.dot(coords_chain1, rotation_matrix.T) + translation
            total_rmsd += np.sum((coords_chain2 - transformed_coords_chain1) ** 2) # Sum of squared distances
        return total_rmsd
    # Initial guess for rotation matrix (identity matrix) and translation vector (zero vector)
    initial_guess = np.zeros(12)
    initial_guess[:9] = np.eye(3).flatten()
    result = minimize(objective, initial_guess, method='BFGS') # Minimize the objective function to estimate rotation and translation
    # Extract rotation and translation from the result
    rotation_matrix = np.reshape(result.x[:9], (3, 3))
    translation = result.x[9:]
    
    return rotation_matrix, translation

def calculate_cross_sectional_area(pdb_file, probe_size=0.6):

    def project_onto_surface(vector_to_project, surface_normal):
        if np.linalg.norm(surface_normal) == 0:
            return vector_to_project
        else:
            surface_normal = surface_normal / np.linalg.norm(surface_normal)
        scalar_projection = np.dot(vector_to_project, surface_normal)
        vector_projection = scalar_projection * surface_normal
        projected_vector = vector_to_project - vector_projection
        return projected_vector

    # Load the GRO file using MDTraj
    traj = mdtraj.load(pdb_file)

    # Select the protein atoms
    protein_traj = traj.atom_slice(traj.top.select('protein'))

    # Calculate the van der Waals surface using MDTraj for the protein
    vdw_surface = mdtraj.shrake_rupley(protein_traj, probe_radius=probe_size)

    R, t = calculate_average_rotation_translation(pdb_file)
    fibril_axis = t
    # Project the van der Waals surface onto the fibril axis
    projected_surface = np.array([project_onto_surface(v, fibril_axis) for v in vdw_surface[0]])
    # Calculate cross-sectional area by summing the projected surface
    cross_section_area = np.sum(projected_surface)

    return cross_section_area


def calculate_average_helical_parameters(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    chain_centroids = {}
    chain_coords = {}

    for model in structure:
        for chain in model:
            coords = []
            for residue in chain:
                for atom in residue:
                    coords.append(atom.get_coord())
            coords = np.array(coords)
            centroid = coords.mean(axis=0)
            chain_centroids[chain.id] = centroid
            chain_coords[chain.id] = coords

    # Infer fibril axis using PCA on centroids
    centroid_matrix = np.array(list(chain_centroids.values()))
    pca = PCA(n_components=1)
    fibril_axis = pca.fit(centroid_matrix).components_[0]
    projections = {cid: np.dot(c, fibril_axis) for cid, c in chain_centroids.items()}

    # Sort chain IDs by projection along the main axis
    sorted_chains = sorted(projections, key=projections.get)

    # Create coordinate pairs of spatial neighbors
    chains_coords = []
    for i in range(len(sorted_chains) - 1):
        c1 = sorted_chains[i]
        c2 = sorted_chains[i + 1]
        coords1 = chain_coords[c1]
        coords2 = chain_coords[c2]
        if coords1.shape[0] == coords2.shape[0]:  # Sanity check
            chains_coords.append((coords1, coords2))

    def objective(params):
        rotation_angle = params[0]
        translation = params[1]
        axis = params[2:5]
        axis = axis / np.linalg.norm(axis)

        # Rodrigues rotation matrix
        ux, uy, uz = axis
        c = np.cos(rotation_angle)
        s = np.sin(rotation_angle)
        C = 1 - c
        R = np.array([
            [c + ux**2*C, ux*uy*C - uz*s, ux*uz*C + uy*s],
            [uy*ux*C + uz*s, c + uy**2*C, uy*uz*C - ux*s],
            [uz*ux*C - uy*s, uz*uy*C + ux*s, c + uz**2*C]
        ])

        total_squared_error = 0.0
        total_points = 0
        for coords1, coords2 in chains_coords:
            transformed = np.dot(coords1, R.T) + translation * axis
            diff = coords2 - transformed
            total_squared_error += np.sum(diff ** 2)
            total_points += coords1.shape[0]
        return total_squared_error / total_points  # MSE

    # Initial guess
    initial_guess = np.zeros(5)
    initial_guess[0] = 0.0  # rotation
    initial_guess[1] = 1.0  # translation guess (e.g., 1 Å)
    initial_guess[2:5] = fibril_axis  # axis guess from PCA

    result = minimize(objective, initial_guess, method='BFGS')

    rotation_angle = result.x[0]
    translation = result.x[1]
    axis = result.x[2:5] / np.linalg.norm(result.x[2:5])
    rmsd = np.sqrt(objective(result.x))

    return rotation_angle, translation, axis, rmsd

def build_fibril(pdb_file, rotation_angle, translation, axis, n_units, output_file):
    """
    Builds a fibril by applying rotations and translations to a chain from a given PDB file
    and assembling the results using Biopython.

    Parameters:
    - pdb_file: Path to the PDB file containing the initial chain.
    - rotation_angle: Rotation angle in radians.
    - translation: Translation distance along the axis of symmetry.
    - axis: 3-element numpy array representing the axis of symmetry (unit vector).
    - n_units: Number of units to assemble in the fibril.
    - output_file: Output PDB file path for the fibril structure.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('input_structure', pdb_file)
    model = structure[0]
    axis = axis / np.linalg.norm(axis)

    # Get all chain centroids and find the one lowest along the axis
    chain_centroids = {}
    for chain in model:
        coords = [atom.get_coord() for residue in chain for atom in residue]
        centroid = np.mean(coords, axis=0)
        chain_centroids[chain.id] = centroid

    projections = {cid: np.dot(centroid, axis) for cid, centroid in chain_centroids.items()}
    sorted_chain_id = sorted(projections, key=projections.get)[0]
    base_chain = model[sorted_chain_id]

    # Create a new structure and model for the output fibril
    fibril_structure = Structure.Structure('fibril_structure')
    fibril_model = Model.Model(0)
    fibril_structure.add(fibril_model)

    def get_rotation_matrix(angle, axis):
        axis = axis / np.linalg.norm(axis)
        c, s = np.cos(angle), np.sin(angle)
        ux, uy, uz = axis
        return np.array([
            [c + ux**2*(1-c), ux*uy*(1-c) - uz*s, ux*uz*(1-c) + uy*s],
            [uy*ux*(1-c) + uz*s, c + uy**2*(1-c), uy*uz*(1-c) - ux*s],
            [uz*ux*(1-c) - uy*s, uz*uy*(1-c) + ux*s, c + uz**2*(1-c)]
        ])

    for i in range(n_units):
        rot_matrix = get_rotation_matrix(rotation_angle * i, axis)
        new_chain = Chain.Chain(chr(65 + i))  # A, B, C, ...

        for residue in base_chain:
            new_residue = Residue.Residue(residue.id, residue.resname, residue.segid)
            for atom in residue:
                new_coord = np.dot(atom.coord, rot_matrix.T) + translation * i * axis
                new_atom = Atom.Atom(
                    atom.name, new_coord, atom.bfactor, atom.occupancy,
                    atom.altloc, atom.fullname, atom.element
                )
                new_residue.add(new_atom)
            new_chain.add(new_residue)

        fibril_model.add(new_chain)

    io = PDBIO()
    io.set_structure(fibril_structure)
    io.save(output_file)


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

import numpy as np
from Bio.PDB import PDBParser, PDBIO

def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2 """
    a = vec1 / np.linalg.norm(vec1)
    b = vec2 / np.linalg.norm(vec2)
    v = np.cross(a, b)
    c = np.dot(a, b)
    if np.isclose(c, -1.0):
        # Vectors are opposite
        orthogonal = np.array([1, 0, 0]) if not np.allclose(a, [1, 0, 0]) else np.array([0, 1, 0])
        v = np.cross(a, orthogonal)
        v = v / np.linalg.norm(v)
        H = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = -np.eye(3) + 2 * np.outer(v, v)
    else:
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]],
                         [v[2], 0, -v[0]],
                         [-v[1], v[0], 0]])
        R = np.eye(3) + kmat + kmat @ kmat * ((1 - c) / (s ** 2))
    return R

def align_axis_to_z(input_pdb, output_pdb, axis):
    # Normalize the growth axis and compute rotation matrix
    axis = axis / np.linalg.norm(axis)
    target_axis = np.array([0, 0, 1])
    rot_matrix = rotation_matrix_from_vectors(axis, target_axis)

    # Parse the PDB
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("fibril", input_pdb)

    # Rotate atoms
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coord = atom.get_coord()
                    new_coord = np.dot(rot_matrix, coord)
                    atom.set_coord(new_coord)

    # Save the rotated structure
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)


def solvate_fibril(
    input_structure,
    output_gro,
    topol_file,
    box=None,
    water_model='spc216.gro',
    gmx_path='gmx',
    expand_box=True
):
    """
    Solvate a built fibril using GROMACS.

    Args:
        input_structure (str): Path to the input GRO or PDB file.
        output_gro (str): Path to the output solvated GRO file.
        topol_file (str): Path to the topology file (topol.top).
        box (tuple or list, optional): Box dimensions (nm) as (x, y, z). If None, box is not changed.
        water_model (str): Water model file for solvation (default: 'spc216.gro').
        gmx_path (str): Command to run GROMACS (e.g., '. ~/load_gromacs.sh; gmx').
        expand_box (bool): Whether to expand the box if box is provided (default: True).
    """
    import os
    structure_for_solvation = input_structure
    if box is not None and expand_box:
        expanded = 'expanded_for_solvation.gro'
        editconf_cmd = f"{gmx_path} editconf -f {input_structure} -o {expanded} -box {box[0]} {box[1]} {box[2]}"
        if os.system(editconf_cmd) != 0:
            raise RuntimeError(f'Failed to run: {editconf_cmd}')
        structure_for_solvation = expanded
    solvate_cmd = f"{gmx_path} solvate -cp {structure_for_solvation} -cs {water_model} -p {topol_file} -o {output_gro}"
    if os.system(solvate_cmd) != 0:
        raise RuntimeError(f'Failed to run: {solvate_cmd}')
    if structure_for_solvation != input_structure:
        os.remove(structure_for_solvation)

def neutralize_system(
    structure_file,
    topol_file,
    mdp_file,
    output_tpr,
    output_gro,
    gmx_load_command=None,
    p_name='NA',
    n_name='CL',
    group_name='SOL',
    maxwarn=1
):
    """
    Prepare the system for ion addition and neutralize using GROMACS.

    Args:
        structure_file (str): Path to the input structure file (usually .gro).
        topol_file (str): Path to the topology file (topol.top).
        mdp_file (str): Path to the ions.mdp file.
        output_tpr (str): Path to the output tpr file.
        output_gro (str): Path to the output neutralized gro file.
        gmx_path (str): Command to run GROMACS (default: 'gmx').
        p_name (str): Name for positive ion (default: 'NA').
        n_name (str): Name for negative ion (default: 'CL').
        group_name (str): Name of the solvent group (default: 'SOL').
        maxwarn (int): Maximum warnings for grompp (default: 1).
    """
    import os
    grompp_cmd = f"{gmx_load_command};gmx grompp -f {mdp_file} -c {structure_file} -p {topol_file} -o {output_tpr} -maxwarn {maxwarn}"
    if os.system(grompp_cmd) != 0:
        raise RuntimeError(f'Failed to run: {grompp_cmd}')
    genion_cmd = f"{gmx_load_command};echo {group_name} | gmx genion -s {output_tpr} -o {output_gro} -p {topol_file} -pname {p_name} -nname {n_name} -neutral"
    if os.system(genion_cmd) != 0:
        raise RuntimeError(f'Failed to run: {genion_cmd}')

def center_fibril(input_file, output_file, gmx_path='gmx'):
    """
    Center the fibril in the simulation box using GROMACS editconf.

    Args:
        input_file (str): Path to the input structure file (PDB or GRO).
        output_file (str): Path to the output centered structure file.
        gmx_path (str): Command to run GROMACS (default: 'gmx').
    """
    import os
    cmd = f"{gmx_path} editconf -f {input_file} -o {output_file} -c"
    if os.system(cmd) != 0:
        raise RuntimeError(f'Failed to run: {cmd}')

def atomtype_protein(input_pdb, output_gro, forcefield, gmx_path='gmx', water_model='spce'):
    """
    Generate a GROMACS-compatible GRO file from a PDB file using pdb2gmx.

    Args:
        input_pdb (str): Path to the input PDB file.
        output_gro (str): Path to the output GRO file.
        forcefield (str): Forcefield name for pdb2gmx.
        gmx_path (str): Command to run GROMACS (default: 'gmx').
        water_model (str): Water model to use (default: 'spce').
    """
    import os
    cmd = f"{gmx_path} pdb2gmx -f {input_pdb} -o {output_gro} -ff {forcefield} -water {water_model} -ignh"
    if os.system(cmd) != 0:
        raise RuntimeError(f'Failed to run: {cmd}')
