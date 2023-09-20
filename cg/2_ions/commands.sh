gmx grompp -f ions.mdp -c ../sa8_cg_20_water.gro -p ../topol.top -o ions.tpr
gmx genion -s ions.tpr -o sa8_cg_20_water_ions.gro -p ../topol.top -pname NA+ -nname CL- -neutral
