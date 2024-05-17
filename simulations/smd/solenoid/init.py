#init_project.py
import signac
import numpy as np
import itertools

#organize based on statepoints which simulations to run

project = signac.init_project('./')
statepoints = [
    {
        "pdbID": "1m8n",
        "region": [19, 98],
        "chain": "A",
        "pull_group1": [19, 41],
        "pull_group2": [70,98],
        "eq_steps": 15000,
        "pull_constant": 5000, # kJ/mol/nm^2
        "pull_rate": 0.01, # nm/ps
        "pull_steps": 250000,
        "data_dir": "/data/yang_lab/nehilpkd/SolenoidBank",
    },
   
]

if __name__ == "__main__":
    if not project.doc.get("response_table"):
        project.doc["response_table"] = {}

    for sp in statepoints:
        if not project.find_jobs(sp):
            print("Initializing job", sp)
            project.open_job(sp).init()
