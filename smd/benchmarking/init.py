#init_project.py
import signac
import numpy as np
import itertools


project = signac.init_project('./')
statepoint = {
    "pdbID": "6Y1A",
    "n_repeats": 5,
    "fiberverse_directory": "/data/yang_lab/nehilpkd/fibrilverse_rcsb",
    "type": "single",
    "pull_chains": [0, 1],
}
# pull_constants = [500, 1000, 2500, 5000, 10000]
# pull_rates = [0.001, 0.01, 0.05, 0.1, .25, .5, 1]
# pull_steps = [250000]
pull_constants = [2500]
pull_rates = [.1]
pull_steps = [250000]

statepoints = []
for pull_constant, pull_rate, pull_steps in itertools.product(pull_constants, pull_rates, pull_steps):
    sp = statepoint.copy()
    sp["pull_constant"] = pull_constant
    sp["pull_rate"] = pull_rate
    sp["pull_steps"] = pull_steps
    statepoints.append(sp)
print(statepoints)

if __name__ == "__main__":
    if not project.doc.get("response_table"):
        project.doc["response_table"] = {}

    for sp in statepoints:
        if not project.find_jobs(sp):
            print("Initializing job", sp)
            project.open_job(sp).init()
        