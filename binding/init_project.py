#init_project.py
import signac
import numpy as np
import itertools

#organize based on statepoints which simulations to run

project = signac.init_project('./')
statepoints = [
    {"sequences":["FGAILSS", "FGAILSS"], "collective_variable":["anti-parallel"], "parallel_bias":False},
    {"sequences":["AAAAAA","AAAAAA"], "collective_variable":["anti-parallel"], "parallel_bias":False},
    {"sequences":["KLVFFAE", "KLVFFAE"], "collective_variable":["anti-parallel"], "parallel_bias":False},
    {"sequences":["GDVIEV", "GDVIEV"], "collective_variable":["anti-parallel"], "parallel_bias":False},
    {"sequences":["AAAAAAAAAAAA","AAAAAAAAAAAA"], "collective_variable":["anti-parallel"], "parallel_bias":False},
]

if __name__ == "__main__":
    if not project.doc.get("response_table"):
        project.doc["response_table"] = {}

    for sp in statepoints:
        print("Initializing job", sp)
        project.open_job(sp).init()