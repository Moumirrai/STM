import argparse
import os

from models import TrussData

from solver import TrussSolver
from plotter import export_vtk
from structure_parser import parse_json_file

# Argument parser
command_parser = argparse.ArgumentParser(description="Solve a truss structure from provided JSON file")

command_parser.add_argument(
    "file_path",
    type=str,
    nargs='?',
    default="./data/star.json",
    help="Path to the input JSON"
)

args = command_parser.parse_args()
input_file_path = args.file_path

if not os.path.exists(input_file_path):
    print(f"Error: Input file not found at '{input_file_path}'")
    exit(1)

print(f"Loading data from: {input_file_path}")

truss: TrussData = parse_json_file(input_file_path)

solver = TrussSolver(truss)

solver.solve()

export_vtk(truss)

exit(0)
