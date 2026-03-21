import numpy as np
from models import TrussData
from solver_debug import TrussSolver
from structure_parser import parse_json_file

np.set_printoptions(
    linewidth=250,
)

truss: TrussData = parse_json_file("./data/tri.json")

solver = TrussSolver(truss)

res = solver.solve()