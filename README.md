BioMatSims
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/biomatsims/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/biomatsims/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/BioMatSims/branch/main/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/BioMatSims/branch/main)



# Usage

## Fibril analysis and extension
```python
from biomatsims.geometry_analysis import calculate_average_rotation_translation
from biomatsims.build import build_fibril
import mbuild

pdb_file = 'some_fibril.pdb'
average_rotation, average_translation = calculate_average_rotation_translation(pdb_file)

mol = mbuild.load(pdb_file)
chain = mol.children[0] # assuming the first child is the chain of interest
predicted_fibril = build_fibril(chain, average_rotation, average_translation, n_chains=40)
```

## Mutate your fibril
```python
from enzy_htp.structure import PDBParser
from enzy_htp.mutation import assign_mutant, mutate_stru
from enzy_htp.mutation.api import sync_mutation_over_chains

n_mutants = 1
pdb_file = "your_fibril.pdb"
parser = PDBParser()

for m in range(n_mutants):
    structure = parser.get_structure(pdb_file)
    mutations = assign_mutant(
        structure,
        pattern="r:1[chain A:all]*1R", # this just performs random mutations across the chain
        chain_sync_list=[tuple(structure.chain_names)], # sync mutations across entire fibril
    )
    mutated_structure = mutate_stru(structure, mutations[0])
    parser.save_structure(
        stru=mutated_structure, 
        outfile=f"{pdb_file.split('.')[0]}_mutant_{m}.pdb"
    )
```

## High-throughput simulation 
Alter the fiberverse database location for your machine in`init.py`, create the list of pdb_ids you want to simulate, specify n_replicates, specify pulling conditions, specify the chains you wish to pull. 

After you have initialized your jobs with `python init.py` perform the following operations to run simulations on your system
```bash
python project.py run -o preprocess_pdb
python project.py run -o create_eq_submission
python project.py run -o run_equilibration
```
After minimization and equilibration have finished run:
```bash
python project.py run -o create_pull_mdp
python project.py run -o create_pull_submission
python project.py run -o preprocess_pull
python project.py run -o run_pull
```
Finally, after the pulling simulation has finished run:
```bash
python project.py run -o run_analysis
```

## Post-simulation analysis
After you have finished your simulations you can calculate mechanical properties of interest:
```python
from biomatsims.geometry_analysis import calculate_cross_sectional_area
from biomatsims.characterization import calculate_variable_over_time
from numpy import max, argmax, array

cross_sectional_area = calculate_cross_sectional_area(fibril_axis, 'protofibril.pdb')
cross_sectional_area = cross_sectional_area * 1e-20 # A^2 to m^2

# Calculate the strain in units of nm/nm
time_length = calculate_variable_over_time('pull_pullx.xvg')
length_over_time = array([l for (t, l) in time_length])
strain = (length_over_time - length_over_time[0]) / length_over_time[0]

# Calculate the stress
time_force = calculate_variable_over_time('pull_pullf.xvg')
force_over_time = array([f for (t, f) in time_force]) * (1e-9) * (1/1000) # kJ/mol/nm to N
stress = force_over_time / cross_sectional_area

# Calculate the ultimate tensile strength
ultimate_tensile_strength = max(stress)

# Calculate the elastic modulus, assuming there is no plastic deformation
max_stress_index = argmax(stress) # assuming the linear elastic region ends at maximum stress
max_stress = stress[max_stress_index]
strain_at_max = strain[max_stress_index]
elastic_modulus = max_stress / strain_at_max
```


A repo for building, characterizing, and analyizing protein fibrils.

### Copyright

Copyright (c) 2024, Kieran Nehil-Puleo


#### Acknowledgements
##### NSF GRFP
##### Yang Lab, Vanderbilt University
##### Cummings Lab, Vanderbilt University
