#init_project.py
import signac
import numpy as np
import itertools

#organize based on statepoints which simulations to run

project = signac.init_project('./')
statepoints = [
    {"amyloids":["FGAILSS", "FGAILSS"], "linkers":["SGRGGLGGQGAGGGAGQGGYGGLGSQGT", "SGRGGLGGQGAGGGAGQGGYGGLGSQGT"], "amyloid_location_index":[12,12], "collective_variable":["anti-parallel"], "parallel_bias":False},
    {"amyloids":["AAAAAA", "AAAAAA"], "linkers":["SGRGGLGGQGAGGGAGQGGYGGLGSQGT", "SGRGGLGGQGAGGGAGQGGYGGLGSQGT"], "amyloid_location_index":[12,12], "collective_variable":["anti-parallel"], "parallel_bias":False},
    {"amyloids":["KLVFFAE", "KLVFFAE"], "linkers":["SGRGGLGGQGAGGGAGQGGYGGLGSQGT", "SGRGGLGGQGAGGGAGQGGYGGLGSQGT"], "amyloid_location_index":[12,12], "collective_variable":["anti-parallel"], "parallel_bias":False},
    {"amyloids":["GDVIEV", "GDVIEV"], "linkers":["SGRGGLGGQGAGGGAGQGGYGGLGSQGT", "SGRGGLGGQGAGGGAGQGGYGGLGSQGT"], "amyloid_location_index":[12,12], "collective_variable":["anti-parallel"], "parallel_bias":False},
]

if __name__ == "__main__":
    if not project.doc.get("response_table"):
        project.doc["response_table"] = {}

    for sp in statepoints:
        print("Initializing job", sp)
        project.open_job(sp).init()