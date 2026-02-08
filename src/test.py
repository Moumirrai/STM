from models import TrussData
from solver_lagrange import LagrangeTrussSolver
from plotter import export_vtk
from structure_parser import parse_json_file
import numpy as np
from scipy.optimize import least_squares
from termcolor import colored

np.set_printoptions(
    linewidth=250,
)


truss: TrussData = parse_json_file('./data/grid.json')

solver = LagrangeTrussSolver(truss)

res = solver.solve()

