import Bio
import Bio.SeqUtils
from Bio.PDB import PDBParser, PDBIO, Select
import pymol2
import numpy as np
import os

def single_to_triple(seq):
	seq = [Bio.SeqUtils.IUPACData.protein_letters_1to3[aa] for aa in seq]
	seq = ''.join(seq)
	t = iter(seq)
	trip_seq = ' '.join((a+b+c).upper() for a,b,c in zip(t,t,t))
	return trip_seq

aa_charge = {
    'A': 0,
    'C': 0,
    'D': -1,
    'E': -1,
    'F': 0,
    'G': 0,
    'H': 1,
    'I': 0,
    'K': 1,
    'L': 0,
    'M': 0,
    'N': 0,
    'P': 0,
    'Q': 0,
    'R': 1,
    'S': 0,
    'T': 0,
    'V': 0,
    'W': 0,
    'Y': 0
}


def search_pattern_in_consecutive_lines(file_path, pattern1, pattern2):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i in range(len(lines) - 1):
            line1 = lines[i]
            line2 = lines[i + 1].strip()
            if pattern1 in line1 and pattern2 in line2:
                return i + 1
    
    return False

def add_line_at_line_number(file_path, line_number, line_content):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    lines.insert(line_number, line_content + '\n')

    with open(file_path, 'w') as file:
        file.writelines(lines)

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
            pass

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

def remove_chains(input_pdb, output_pdb, chains_to_remove):
    class ChainSelector(Select):
        def accept_chain(self, chain):
            return chain.id not in chains_to_remove

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=ChainSelector())


def extract_chain(input_pdb, chain_id):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)

    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                return chain

    return None


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




def clean_structure(input_pdb, output_pdb, working_dir=None):
    """
    Clean a PDB structure by removing water and non-standard residues.

    Args:
        input_pdb (str): Path to the input PDB file.
        output_pdb (str): Path to the output cleaned PDB file.
    """
    pdbID = os.path.basename(input_pdb).split('.')[0]

    renamed_pdb = f"{working_dir}/{pdbID}_renamed.pdb"
    renumbered_pdb = f"{working_dir}/{pdbID}_renumbered.pdb"
    cleaned_pdb = f"{working_dir}/{pdbID}_cleaned.pdb"
    rename_amino_acid(input_pdb, renamed_pdb)
    renumber_residues(renamed_pdb, renumbered_pdb)
    remove_ligands_water_ions(renumbered_pdb, cleaned_pdb)
    add_hydrogens(cleaned_pdb, output_pdb)





if __name__ == "__main__":
	res = single_to_triple('GAG')
	print(res)
	 

