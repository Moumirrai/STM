import numpy as np
from vornoi_structure import generate_voronoi_structure
from parameter_solver import solveParameters_iso, solveParameters_orto
from models import TrussData
from plotter import plot_deformed_structure
from solver import TrussSolver
from structure_parser import parse_structure_data, truss_data_to_json_string


np.set_printoptions(
    linewidth=250,
)

trussData = generate_voronoi_structure(10, 5, 1000, 0.2)

solveParameters_iso(trussData)

""" truss: TrussData = parse_structure_data(
    trussData, explicitEigenStrain=np.array([0.0, 0.0, 0.0])
)

#truss: TrussData = parse_json_file("./data/grid.json")

solver = TrussSolver(truss)

res = solver.solve()

plot_deformed_structure(truss,tiled=False)
 """