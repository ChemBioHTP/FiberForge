from Bio.PDB import PDBParser
from scipy.optimize import minimize
import numpy as np
from sklearn.decomposition import PCA
import mdtraj


def calculate_cross_sectional_area(fibril_axis, pdb_file, probe_size=0.6, n_cross_sections=10):

    def project_onto_surface(vector_to_project, surface_normal):
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

    # Identify the indices of atoms that are part of the fibril
    fibril_atoms_indices = np.array([atom.index for atom in traj.topology.atoms if atom.residue.is_protein])

    # Calculate the cross-sectional area for each frame along the fibril axis
    cross_section_areas = []
    for frame_idx in range(traj.n_frames):
        # fibril_frame = traj.xyz[frame_idx][fibril_atoms_indices]
        # axis_frame = traj.xyz[frame_idx][fibril_axis]

        # Project the van der Waals surface onto the fibril axis
        projected_surface = np.array([project_onto_surface(v, fibril_axis) for v in vdw_surface[frame_idx]])

        # Calculate cross-sectional area by summing the projected surface
        cross_section_area = np.sum(projected_surface)

        cross_section_areas.append(cross_section_area)

    # Calculate the average cross-sectional area
    average_cross_section_area = np.mean(cross_section_areas)

    return average_cross_section_area

def estimate_rotation_translation_between_chains(pdb_file, chain1_id, chain2_id):
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
    
    # Define the function to minimize (sum of squared distances)
    def objective(params):
        rotation_matrix = np.reshape(params[:9], (3, 3))
        translation = params[9:]
        transformed_coords_chain1 = np.dot(coords_chain1, rotation_matrix.T) + translation
        return np.sum((coords_chain2 - transformed_coords_chain1) ** 2)
    
    # Initial guess for rotation matrix (identity matrix) and translation vector (zero vector)
    initial_guess = np.zeros(12)
    initial_guess[:9] = np.eye(3).flatten()
    
    # Minimize the objective function to estimate rotation and translation
    result = minimize(objective, initial_guess, method='BFGS', tol=1e-6)
    
    # Extract rotation and translation from the result
    rotation_matrix = np.reshape(result.x[:9], (3, 3))
    translation = result.x[9:]
    
    return rotation_matrix, translation

def calculate_fibril_dimensions(coords):
    pca = PCA(n_components=3)
    pca.fit(coords)
    return pca.explained_variance_

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