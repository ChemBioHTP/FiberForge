BioMatSims
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/biomatsims/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/biomatsims/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/BioMatSims/branch/main/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/BioMatSims/branch/main)



# Example usage

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


A repo for building, characterizing, and analyizing protein fibrils.

### Copyright

Copyright (c) 2024, Kieran Nehil-Puleo


#### Acknowledgements
##### NSF GRFP
##### Yang Lab, Vanderbilt University
##### Cummings Lab, Vanderbilt University
