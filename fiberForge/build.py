from Bio.PDB import PDBParser
import numpy as np
from plotly import graph_objects as go
from fiberForge.geometry_analysis import estimate_rotation_translation_between_chains
from scipy.optimize import minimize
import mdtraj


def visualize_rotation_translation(pdb_file, chain1_id, chain2_id):
    rotation, translation = estimate_rotation_translation_between_chains(pdb_file, chain1_id, chain2_id)
    
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    
    fig = go.Figure()
    
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
    
    # Plot the original and transformed coordinates
    x_chain2, y_chain2, z_chain2 = coords_chain2[:,0], coords_chain2[:,1], coords_chain2[:,2]
    x_transformed_chain1, y_transformed_chain1, z_transformed_chain1 = transformed_coords_chain1[:,0], transformed_coords_chain1[:,1], transformed_coords_chain1[:,2]
    
    # fig.add_trace(go.Scatter3d(x=x_chain1, y=y_chain1, z=z_chain1, mode='markers', name='Chain 1 (Original)'))
    fig.add_trace(go.Scatter3d(x=x_chain2, y=y_chain2, z=z_chain2, mode='markers', name='Chain 2'))
    fig.add_trace(go.Scatter3d(x=x_transformed_chain1, y=y_transformed_chain1, z=z_transformed_chain1, mode='markers', name='Chain 1 (Transformed)'))
    
    fig.update_layout(title='Rotation and Translation Visualization',
                      scene=dict(
                          xaxis_title='X',
                          yaxis_title='Y',
                          zaxis_title='Z'
                      ))
    
    fig.show()

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

def calculate_cross_sectional_area(job, probe_size=0.6):

    def project_onto_surface(vector_to_project, surface_normal):
        surface_normal = surface_normal / np.linalg.norm(surface_normal)
        scalar_projection = np.dot(vector_to_project, surface_normal)
        vector_projection = scalar_projection * surface_normal
        projected_vector = vector_to_project - vector_projection
        return projected_vector
    
    pdb_file = job.path + '/0_preprocess/protofibril.pdb'

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

def build_fibril(chain, rotation, translation, n_units):
    import mbuild as mb
    predicted_fibril = mb.Compound()
    for i in range(n_units):
        chain_copy = mb.clone(chain)
        rotation_matrix = np.linalg.matrix_power(rotation, i)
        chain_copy.xyz = np.dot(chain_copy.xyz, rotation_matrix.T)
        chain_copy.translate(translation * i / 10.0) # need to scale translation and convert to correct units
        predicted_fibril.add(chain_copy)
    return predicted_fibril

def build_fibril_biopython(pdb_file, rotation, translation, n_units, output_file):
    """
    Builds a fibril by applying rotations and translations to a chain from a given PDB file
    and assembling the results using Biopython.

    Parameters:
    - pdb_file: The path to the PDB file containing the initial chain.
    - rotation: A 3x3 numpy array representing the rotation matrix.
    - translation: A 3-element numpy array representing the translation vector.
    - n_units: The number of units to assemble in the fibril.
    - output_file: The path to the output PDB file where the fibril structure will be saved.
    
    Returns:
    - None: The fibril structure is saved to the output PDB file.
    """
    import numpy as np
    from Bio.PDB import PDBParser, PDBIO, Superimposer
    from Bio.PDB import Structure, Model, Chain, Residue, Atom
    
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('input_structure', pdb_file)
    
    # Get the first model, chain, and residue as the starting point
    model = structure[0]
    chain = model.child_list[0]  # Assumes single chain in input
    
    # Create a new structure for the fibril
    fibril_structure = Structure.Structure('fibril_structure')
    fibril_model = Model.Model(0)
    fibril_structure.add(fibril_model)
    
    # Loop over the number of units
    for i in range(n_units):
        # Calculate the rotation matrix for the current unit
        rotation_matrix = np.linalg.matrix_power(rotation, i)
        
        # Create a new chain for the current unit
        new_chain = Chain.Chain(chr(65 + i))  # Chain names: A, B, C, ...
        
        # Copy residues from the original chain and apply transformations
        for residue in chain:
            new_residue = Residue.Residue(residue.id, residue.resname, residue.segid)
            for atom in residue:
                new_atom = Atom.Atom(atom.name, np.dot(atom.coord, rotation_matrix.T) + translation * i, atom.bfactor, atom.occupancy, atom.altloc, atom.fullname, atom.element)
                new_residue.add(new_atom)
            new_chain.add(new_residue)
        
        # Add the new chain to the fibril model
        fibril_model.add(new_chain)
    
    # Save the fibril structure to a PDB file
    io = PDBIO()
    io.set_structure(fibril_structure)
    io.save(output_file)
