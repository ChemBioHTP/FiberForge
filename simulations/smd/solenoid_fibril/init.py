#init_project.py
import signac
import numpy as np
import itertools
from copy import deepcopy

#organize based on statepoints which simulations to run

project = signac.init_project('./')
sp = {
        # "chain_COM_region": "middle",
        # "chain_COM_res_interval": [-5,5],
        "pull_group1": 0,
        "pull_group2": -1,
        "eq_steps": 15000,
        "pull_constant": 5000, # kJ/mol/nm^2
        "pull_rate": 0.01, # nm/ps
        "pull_steps": 250000,
        "data_dir": "/data/yang_lab/nehilpkd/fiberForge/smd/solenoid_fibril/dataset",
}

pdbs = [
    '2rnm',
    '5aej',
    '5wor',
    '6eka',
    '6ria',
    '6rib',
    '7b90',
    '7bfc',
    '7pnb',
    '7pqc',
    '7y8q',
    '8cio'
]
statepoints = []

if __name__ == "__main__":
    if not project.doc.get("response_table"):
        project.doc["response_table"] = {}

    for pdb in pdbs:
        new_sp = deepcopy(sp)
        new_sp['pdbID'] = pdb
        statepoints.append(new_sp)

    for sp in statepoints:
        if not project.find_jobs(sp):
            print("Initializing job", sp)
            project.open_job(sp).init()
