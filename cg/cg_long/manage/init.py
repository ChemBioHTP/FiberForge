import signac

project = signac.init_project()

for p in range(1, 10):
    '''
    data space := (protein_seq, number of proteins, )
    '''
    sp = {"p": p, "kT": 1.0, "N": 1000}
    job = project.open_job(sp).init()