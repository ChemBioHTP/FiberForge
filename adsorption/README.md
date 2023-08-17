# Procedure
following the tutorial here https://ambermd.org/tutorials/advanced/tutorial27/pro_metal.php

1. First step convert .car file to .pdb with VMD so the CRYST1 information is correct, the raw .pdb file from InterfaceMD doesn't contain this information 
2. Next we tile the silica surface to the desired shape with PropPDB from Amber
3. Next we copy over the mfp pdb from alphafold and process it
    1.  the first line (MODEL) 
    2. Remove the second to last line (ENDMDL)
Note I manually had to move with VMD the mfp and silica surface to positions where they don't intersect since this causes vmd to infer bonding and makes you unable to move just one of the molecules
4. use cat to combine the surface and the protein
4. Finally we  use vmd 

