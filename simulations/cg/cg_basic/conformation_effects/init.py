import signac
from pymol import cmd

project = signac.init_project()

# mol = pymol.cmd.load('')

# def make_mutation(p, c, r):
#     rotkit.mutate(p, chain=c, resi=r, target="CYS", mutframe=1)
#     cmd.select("%s%s%s"%(p,c,r),"/%s//%s/%s"%((p,c,r)))
#     cmd.save('mut.pdb')

conformations = list(range(10))

for conformer in range(len(conformations)):
    '''
    data space := (protein_seq, number of proteins, )
    '''
    sp = {"conformer": conformer}
    job = project.open_job(sp).init()