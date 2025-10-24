FiberForge
==============================

# Installation

```python
conda install -f environment.yml
conda activate ff
pip install .
```


# Usage

## Fibril analysis and extension
```python
from FiberForge.build import (
  identify_protofibrils,
  calculate_average_helical_parameters,
  build_fibril,
  solvate_fibril
)
from FiberForge.utils import (
    clean_structure,
    remove_chains
)

# The path to your amyloid structure in PDB format
pdb_file = "amyloid.pdb" # A classic is 2MXU (ABeta42)

# Locate the protofibrils, returns List[dict[str, np.array(3)]]
protofibrils = identify_protofibrils(pdb_file)

# If multiple protofibril structure, remove other chains                    
remove_chains(
    pdb_file, 
    "protofibril.pdb", # Saves to protofibril.pdb
    chains_to_remove = [] # for example we assume fibril made of 1 protofibril
) 

# Calculate the best theta and t 
theta, t, growth_axis = calculate_average_helical_parameters("protofibril.pdb")

build_fibril(
    "protofibril.pdb", 
    theta, 
    t, 
    n_units=20, 
    output_file = "extended_protofibril.pdb"
)

solvate_fibril(
    input_structure="extended_protofibril.gro",
    output_gro="solvated.gro",
    topol_file="topol.top",
    box=box,
    water_model='spc216.gro',
    gmx_path='module load gromacs;gmx',
    expand_box=False
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
from FiberForge.analyze import calculate_cross_sectional_area

pdb_file = 'extended_protofibril.pdb'
project_path = 'path/to/files'

cross_sectional_area = calculate_cross_sectional_area(pdb_file)
cross_sectional_area = cross_sectional_area * 1e-20 # A^2 to m^2

# Calculate the strain in units of nm
time_length = calculate_variable_over_time(project_path + '/4_smd/pull_pullx.xvg')
length_over_time = np.array([l for (t, l) in time_length])
strain = (length_over_time - length_over_time[0]) / length_over_time[0]

# Calculate the stress
time_force = calculate_variable_over_time(project_path + '/4_smd/pull_pullf.xvg')
force_over_time = np.array([f for (t, f) in time_force])*(1e-9)*(1/1000) # kJ/mol/nm to N
stress = force_over_time / cross_sectional_area

# Calculate the ultimate tensile strength
ultimate_tensile_strength = np.max(stress)

# Calculate the elastic modulus, assuming there is no plastic deformation
E, yield_point = estimate_elastic_modulus(stress, strain)
```


A repo for building, characterizing, and analyizing protein fibrils.

### Copyright

Copyright (c) 2025, Kieran Nehil-Puleo


#### Acknowledgements
##### NSF GRFP
##### Yang Lab, Vanderbilt University
##### Cummings Lab, Vanderbilt University
