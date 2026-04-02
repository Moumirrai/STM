import numpy as np
from joblib import Parallel, delayed
from vornoi_structure import generate_voronoi_structure
from parameter_solver import solveParameters_iso

num_tries = 80
x = np.logspace(0, -1, 20)  # logaritmické rozložení

def single_run(x_i):
    trussData = generate_voronoi_structure(5, 5, 500, x_i)
    return solveParameters_iso(trussData)

datapoints = []

for x_i in x:
    # n_jobs=-1 použije všechna CPU jádra
    results = Parallel(n_jobs=-1)(
        delayed(single_run)(x_i) for _ in range(num_tries)
    )
    E, v = zip(*results)
    datapoints.append((x_i, np.mean(E), np.std(E), np.mean(v), np.std(v)))
    print(f"Hotovo x={x_i:.3f}")

np.savetxt("vornoi_fidelity.csv", np.array(datapoints), delimiter=",", comments="")