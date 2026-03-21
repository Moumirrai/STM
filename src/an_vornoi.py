import numpy as np
from models import TrussData
from solver_debug import TrussSolver
from structure_parser import parse_json_file
from plotter import plot_deformed_structure
from parameter_solver import solveParameters_iso

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

solveParameters_iso(f"./structures/v{json_number}.json")