import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue, Atom
import signac
from scipy.optimize import minimize
from fiberForge.dataset import pdbs
from fiberForge.build import identify_protofibrils, build_fibril_from_scalar_with_biopython, calculate_average_helical_parameters
from fiberForge.utils import *
from Bio.PDB import PDBParser, Superimposer
import pickle
import warnings
from tqdm import tqdm
warnings.filterwarnings("ignore")
import scienceplots
import matplotlib.pyplot as plt
plt.style.use([
    "science",
    'no-latex'
])

rmsds = {}
rotations = {}
translations = {}
for i, p in tqdm(enumerate(pdbs), total=len(pdbs)):
    try:
        pdb_path = '/panfs/accrepfs.vampire/data/yang_lab/nehilpkd/fibrilverse_rcsb' + '/'+ p + '.pdb'
        os.system(f"cp {pdb_path} /panfs/accrepfs.vampire/data/yang_lab/nehilpkd/fiberForge/notebooks/tmp3_path_for_pdbs")
        pdb_path = "/panfs/accrepfs.vampire/data/yang_lab/nehilpkd/fiberForge/notebooks/tmp3_path_for_pdbs/" + p + ".pdb"
        cleaned_pdb = f"/panfs/accrepfs.vampire/data/yang_lab/nehilpkd/fiberForge/notebooks/tmp3_path_for_pdbs/cleaned_{p}.pdb"
        # Clean the chains
        remove_ligands_water_ions(pdb_path, cleaned_pdb)

        protofibrils = identify_protofibrils(cleaned_pdb, 8.0) # 8.0 is the distance cutoff for identifying protofibrils
        if len(protofibrils) > 1:
            chains_to_remove = []
            for i in range(1, len(protofibrils)):
                chains_to_remove.extend([chain_id for chain_id in protofibrils[i].keys()])
            remove_chains(
                cleaned_pdb, 
                f"/data/yang_lab/nehilpkd/fiberForge/notebooks/tmp3_path_for_pdbs/exp_{p}_protofibril.pdb", 
                chains_to_remove
            )
        else:
            os.rename(cleaned_pdb, f"/panfs/accrepfs.vampire/data/yang_lab/nehilpkd/fiberForge/notebooks/tmp3_path_for_pdbs/exp_{p}_protofibril.pdb")
        pdb_file = f"/panfs/accrepfs.vampire/data/yang_lab/nehilpkd/fiberForge/notebooks/tmp3_path_for_pdbs/exp_{p}_protofibril.pdb"
        rotation_angle, translation, axis, rmsd = calculate_average_helical_parameters(pdb_file)
        rotations[p] = rotation_angle
        translations[p] = translation
        build_fibril_from_scalar_with_biopython(
            pdb_file = pdb_file, 
            rotation_angle = rotation_angle, 
            translation = translation,
            axis = axis,
            n_units = len(protofibrils[0]),
            output_file = f"/panfs/accrepfs.vampire/data/yang_lab/nehilpkd/fiberForge/notebooks/tmp3_path_for_pdbs/recon_{p}_protofibril.pdb"
        )

        def calculate_total_rmsd(structure1, structure2):
            parser = PDBParser(QUIET=True)
            s1 = parser.get_structure('model1', structure1)
            s2 = parser.get_structure('model2', structure2)
            
            atoms1 = []
            atoms2 = []
            num_chains = 0

            for chain1, chain2 in zip(s1[0], s2[0]):  # Assuming only one model in each PDB
                chain_atoms1 = list(chain1.get_atoms())
                chain_atoms2 = list(chain2.get_atoms())
                
                if len(chain_atoms1) != len(chain_atoms2):
                    print(f"Chain {chain1.id} length mismatch, skipping RMSD calculation for this chain.")
                    continue
                
                atoms1.extend(chain_atoms1)
                atoms2.extend(chain_atoms2)
                num_chains += 1
            
            if num_chains == 0:
                raise ValueError("No valid chains found for RMSD calculation.")
            
            # Perform superimposition
            super_imposer = Superimposer()
            super_imposer.set_atoms(atoms1, atoms2)
            super_imposer.apply(s2[0].get_atoms())  # Apply transformation to structure 2
            
            # Compute RMSD
            total_rmsd = super_imposer.rms
            per_chain_rmsd = total_rmsd / num_chains

            return total_rmsd, per_chain_rmsd

        # Usage
        total_rmsd, per_chain_rmsd = calculate_total_rmsd(
            pdb_file,
            f"/panfs/accrepfs.vampire/data/yang_lab/nehilpkd/fiberForge/notebooks/tmp3_path_for_pdbs/recon_{p}_protofibril.pdb"
        )

        rmsds[p] = per_chain_rmsd

    except Exception as e:
        print(f"Error processing job {p}: {e}")
        continue


# pickle the rmsds
with open('/panfs/accrepfs.vampire/data/yang_lab/nehilpkd/fiberForge/notebooks/tmp3_path_for_pdbs/rmsds.pkl', 'wb') as f:
    pickle.dump(rmsds, f)
with open('/panfs/accrepfs.vampire/data/yang_lab/nehilpkd/fiberForge/notebooks/tmp3_path_for_pdbs/rotations.pkl', 'wb') as f:
    pickle.dump(rotations, f)
with open('/panfs/accrepfs.vampire/data/yang_lab/nehilpkd/fiberForge/notebooks/tmp3_path_for_pdbs/translations.pkl', 'wb') as f:
    pickle.dump(translations, f)