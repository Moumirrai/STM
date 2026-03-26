import numpy as np
import re
from pathlib import Path
from vornoi_structure import generate_voronoi_structure
from parameter_solver import solveParameters_iso, solveParameters_orto
from models import TrussData
from plotter import plot_deformed_structure
from solver import TrussSolver
from structure_parser import parse_structure_data, truss_data_to_json_string

num_tries = 80
x = np.linspace(1, 0.1, 20, endpoint=True)

datapoints = []

for x_i in x:
    iter_res = []
    for inner_i in range(num_tries):
        trussData = generate_voronoi_structure(5, 5, 500, x_i)
        iter_res.append(solveParameters_iso(trussData))
        print(f"Completed {inner_i+1}/{num_tries} for x={x_i:.3f}", end="\r", flush=True)
        
    #zpočítej střední hodnotu a smerodatnou odchylku pro každý parametr, máme pouze E a v
    E, v = zip(*iter_res)
    E_mean = np.mean(E)
    E_std = np.std(E)
    v_mean = np.mean(v)
    v_std = np.std(v)
    
    datapoints.append((x_i, E_mean, E_std, v_mean, v_std))
    
#save to csv
np.savetxt(
    "vornoi_fidelity.csv", np.array(datapoints), delimiter=",", comments="",
)
    