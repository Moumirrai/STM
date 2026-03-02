import matplotlib
import math

import numpy as np

from generator import create_tie_structure, create_tie_structure_angle
from models import TrussData
from parameter_solver import solveParameters_iso, solveParameters_orto
from plotter import plot_deformed_structure

#matplotlib.use("QtAgg")
import matplotlib.pyplot as plt

from solver import TrussSolver
from structure_parser import parse_structure_data

height = 0.1
width = height * 1
max_angle = math.degrees(math.atan(height / width))

angle = 20

structure = create_tie_structure_angle(height, width, angle)

solveParameters_orto(create_tie_structure_angle(height, width, angle))

truss: TrussData = parse_structure_data(
    structure, explicitEigenStrain=np.array([1, 0, 0]),
)
solver = TrussSolver(truss)
res = solver.solve()
plot_deformed_structure(truss)


