#init_project.py
import signac
import numpy as np
import itertools

#organize based on statepoints which simulations to run

project = signac.init_project('./')
statepoints = [
    {"sequence":"FGAILS", 'n_molecules': 30},
    {"sequence":"VKVKVKVKVPPTKVKVKVKV", 'n_molecules': 10},
    {"sequence": "VKVKVEVK", 'n_molecules': 30},
    {"sequence": "FEFEFKFK", 'n_molecules': 30},
    {"sequence": "ADARADARADARADA", 'n_molecules': 12},
    {"sequence": "AEAEAKAKAEAEAKAK", 'n_molecules': 30},

   
]

if __name__ == "__main__":
    if not project.doc.get("response_table"):
        project.doc["response_table"] = {}
    
    already_initalized_jobs = [j.sp for j in project]
    for sp in statepoints:
        if sp not in already_initalized_jobs:
            print("Initializing job", sp)
            project.open_job(sp).init()