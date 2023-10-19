# init.py
import signac

project = signac.init_project()

for p in range(1, 10):
    sp = {"p": p, "kT": 1.0, "N": 1000}
    job = project.open_job(sp).init()
