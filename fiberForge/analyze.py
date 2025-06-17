import pandas as pd
import numpy as np
from scipy.interpolate import CubicSpline
from Bio.PDB import PDBParser
import math

def estimate_elastic_modulus(stress, strain, window_size=400,):
    """
    Estimate the elastic modulus from a stress-strain curve.
    
    Parameters
    ----------
    stress : array-like
        The stress values.
    strain : array-like
        The strain values.
    
    Returns
    -------
    float
        The elastic modulus.
    """
    df = pd.DataFrame({'strain': strain, 'stress': stress})
    while True:
        moving_avg = df.rolling(window=window_size).mean()
        if moving_avg.isnull().values.all():
            window_size = window_size // 2
        else:
            break
    moving_avg = moving_avg.dropna()
    

    args = np.argsort(moving_avg['strain']).values
    x = moving_avg['strain'].values[args]
    y = moving_avg['stress'].values[args]

    # Make x strictly increasing
    x, idx = np.unique(x, return_index=True)
    y = y[idx]

    # Create a cubic spline
    cs = CubicSpline(x, y)

    # Calculate the first derivative of the spline
    cs_derivative = cs.derivative()

    # Create a fine grid of x values for smooth plotting
    x_fine = np.linspace(x.min(), x.max(), 2000)
    y_fine = cs(x_fine)

    dy2_dx2_fine = cs_derivative.derivative()(x_fine)

    inflection_points = np.isclose(dy2_dx2_fine, 0.0, atol=1e13) == 1

    E = (y_fine[inflection_points][0] - y_fine[0]) / (x_fine[inflection_points][0] - x_fine[0])

    yield_point = x_fine[inflection_points][0]

    return E, yield_point

def calculate_variable_over_time(xvg_file):
    with open(xvg_file, 'r') as f:
        lines = f.readlines()
        time_data = []
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                time_data.append(list(map(float, line.split())))
    return time_data

def count_interchain_hydrogen_bonds(pdb_file, distance_cutoff=3.5):
    """
    Count the number of hydrogen bonds between chains in an amyloid structure.
    Uses a simple distance-based criterion between N and O atoms of different chains.

    Parameters
    ----------
    pdb_file : str
        Path to the PDB file.
    distance_cutoff : float
        Maximum distance (in Å) between donor and acceptor atoms to be considered a hydrogen bond.

    Returns
    -------
    int
        Number of interchain hydrogen bonds.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('amyloid', pdb_file)
    atoms_by_chain = {}
    for model in structure:
        for chain in model:
            atoms_by_chain[chain.id] = [atom for residue in chain for atom in residue if atom.element in ('N', 'O')]
    chains = list(atoms_by_chain.keys())
    hbonds = 0
    for i, chain1 in enumerate(chains):
        for chain2 in chains[i+1:]:
            for atom1 in atoms_by_chain[chain1]:
                for atom2 in atoms_by_chain[chain2]:
                    if atom1.element != atom2.element:  # N-O or O-N
                        dist = atom1 - atom2
                        if dist <= distance_cutoff:
                            hbonds += 1
    return hbonds

def calculate_interchain_lj_energy(pdb_file, epsilon=0.1, sigma=3.5, cutoff=10.0):
    """
    Calculate the total Lennard-Jones potential energy between chains in an amyloid structure.
    Only considers atom pairs between different chains within a cutoff distance.

    Parameters
    ----------
    pdb_file : str
        Path to the PDB file.
    epsilon : float
        Depth of the potential well (kcal/mol or kJ/mol, depending on units).
    sigma : float
        Finite distance at which the inter-particle potential is zero (Å).
    cutoff : float
        Maximum distance (Å) to consider for LJ interactions.

    Returns
    -------
    float
        Total Lennard-Jones energy between chains.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('amyloid', pdb_file)
    atoms_by_chain = {}
    for model in structure:
        for chain in model:
            atoms_by_chain[chain.id] = [atom for residue in chain for atom in residue if atom.element != 'H']
    chains = list(atoms_by_chain.keys())
    total_energy = 0.0
    for i, chain1 in enumerate(chains):
        for chain2 in chains[i+1:]:
            for atom1 in atoms_by_chain[chain1]:
                coord1 = atom1.coord
                for atom2 in atoms_by_chain[chain2]:
                    coord2 = atom2.coord
                    r = np.linalg.norm(coord1 - coord2)
                    if r <= cutoff and r > 0:
                        lj = 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)
                        total_energy += lj
    return total_energy


