#init_project.py
import signac
import numpy as np
import itertools


project = signac.init_project('./')
statepoint = {
    "pdbID": "6Y1A",
    "n_repeats": 1,
    "fiberverse_directory": "/data/yang_lab/nehilpkd/fibrilverse_rcsb",
    "type": "single",
    "protofibril_distance_threshold" : 20.0,
    "pull_chains": [0, 1],
    "build_fiber": True,
}
# pull_constants = [500, 1000, 2500, 5000, 10000]
# pull_rates = [0.001, 0.01, 0.05, 0.1, .25, .5, 1]
# pull_steps = [250000]
pull_constants = [1000]
pull_rates = [.01]
pull_steps = [250000]
box_sizes = [1, 1.5, 2, 2.5, 3]
# fiber_repeats = [10, 20, 30, 40]
n_replicates = 3

statepoints = []

for pull_constant, pull_rate, pull_steps, box_factor in itertools.product(pull_constants, pull_rates, pull_steps, box_sizes):
    for i in range(n_replicates):
        sp = statepoint.copy()
        sp["n_repeats"] = i
        sp["pull_constant"] = pull_constant
        sp["pull_rate"] = pull_rate
        sp["pull_steps"] = pull_steps
        sp["box_factor"] = box_factor
        statepoints.append(sp)

if __name__ == "__main__":
    if not project.doc.get("response_table"):
        project.doc["response_table"] = {}

    for sp in statepoints:
        if not project.find_jobs(sp):
            print("Initializing job", sp)
            project.open_job(sp).init()
        
