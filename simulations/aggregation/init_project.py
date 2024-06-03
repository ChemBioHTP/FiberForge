#init_project.py
import signac
import numpy as np
import itertools

#organize based on statepoints which simulations to run

project = signac.init_project('./')
statepoints = [
    {"sequence":"FGAILS", 'n_molecules': 30},
    {"sequence":"VKVKVKVKVPPTKVKVKVKV", 'n_molecules': 30},
   
]

if __name__ == "__main__":
    if not project.doc.get("response_table"):
        project.doc["response_table"] = {}

    for sp in statepoints:
        print("Initializing job", sp)
        project.open_job(sp).init()