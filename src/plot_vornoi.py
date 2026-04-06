import numpy as np
from models import TrussData
from solver_debug import TrussSolver
from structure_parser import parse_json_file
from plotter import plot_deformed_structure, export_vtk, plot_deformed_structure_black

np.set_printoptions(
    linewidth=250,
)

#get command argument for json number
import sys
if len(sys.argv) > 1:
    json_number = sys.argv[1]
else:
    print("No JSON number provided. Usage: python debug.py <json_number>")
    sys.exit(1)

truss: TrussData = parse_json_file(f"./structures/v{json_number}.json", explicitEigenStrain=np.array([0.0, 0.0, 0.0]))

solver = TrussSolver(truss)

res = solver.solve()
plot_deformed_structure_black(truss, tiled=True)
export_vtk(truss)