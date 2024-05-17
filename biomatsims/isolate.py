from Bio.PDB import PDBParser
import numpy as np
from sklearn.decomposition import PCA


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