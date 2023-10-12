#num of proteins, num of conformers, 

python martinize.py -f sa8.pdb -o topol.top -x sa8_cg.pdb -ff martini22
gmx insert-molecules -box 15 15 15 -nmol 20 -ci sa8_cg.pdb -radius 0.4 -o sa8_cg_20.gro -try 1000
# if you need to regenerate the equilibrated water box, use command below
insane -x 10 -y 10 -z 10 -d 0 -pbc cubic -sol W -excl -0.5 -o waterbox_10nm.gro
gmx solvate -cp sa8_cg_20.gro -cs eq.gro -radius 0.21 -o sa8_cg_20_water.gro -p topol.top
# for some reason this isn't quite working, need to manually fix topol.top with correct martini.itp and space between molecules listed
