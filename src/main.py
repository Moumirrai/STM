import json
import argparse
import os
from models import Node, Element, TrussData

from plotter import render_truss_structure, visualise_axial_forces
from solver import TrussSolver
from src.plotter import export_vtk


def load_data(file_path: str) -> TrussData:
    with open(file_path) as f:
        data = json.load(f)

    totlal_constraints = 0

    nodes = []
    for node_data in data["nodes"]:
        constraints = node_data.get("constraints", "")
        constrained_x = "x" in constraints
        constrained_y = "y" in constraints

        totlal_constraints += int(constrained_x) + int(constrained_y)

        deformations = node_data.get("deformations", {})
        deformation_x = deformations.get("x", 0.0)
        deformation_y = deformations.get("y", 0.0)

        loads = node_data.get("loads", {})
        load_x = loads.get("x", 0.0)
        load_y = loads.get("y", 0.0)

        node = Node(
            dx=eval(node_data["dx"]),
            dy=eval(node_data["dy"]),
            constrained_x=constrained_x,
            constrained_y=constrained_y,
            deformation_x=deformation_x,
            deformation_y=deformation_y,
            load_x=load_x,
            load_y=load_y,
        )

        nodes.append(node)

    elements = [
        Element(
            nodes=(
                nodes[int(element["starting_node"])],
                nodes[int(element["ending_node"])],
            ),
        )
        for element in data["elements"]
    ]
    print(f"Total constraints: {totlal_constraints}")
    return TrussData(nodes, elements, totlal_constraints)


# Argument parser
parser = argparse.ArgumentParser(description="Solve a truss structure from privided JSON file")

parser.add_argument(
    "file_path", 
    type=str, 
    nargs='?',
    default="../data/default.json",
    help="Path to the input JSON"
)

args = parser.parse_args()
input_file_path = args.file_path

if not os.path.exists(input_file_path):
    print(f"Error: Input file not found at '{input_file_path}'")
    exit(1)

print(f"Loading data from: {input_file_path}")
truss: TrussData = load_data(input_file_path)


#render_truss_structure(truss)

solver = TrussSolver(truss)

solver.solve()

#visualise_axial_forces(truss)
export_vtk(truss)

exit(0)
