import math

import matplotlib

# matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
import numpy as np

from models import TrussData
from plotter import export_vtk, plot_deformed_structure
from solver import TrussSolver
from structure_parser import parse_json_file, parse_structure_data

truss: TrussData = parse_json_file("./data/simple_truss1.json")

solver = TrussSolver(truss)

res = solver.solve()

for element in truss.elements:
    print("Cos sin")
    print(element.get_cos_sin())
    print(element._stiffness_matrix)


exit(0)

export_vtk(truss)

plot_deformed_structure(truss, tiled=False, original=True)
