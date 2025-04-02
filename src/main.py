import json
import matplotlib.pyplot as plt
import numpy as np

from models import Node, Element, TrussData

from plotter import render_truss_structure
from solver import TrussSolver


def load_data(file_path: str) -> TrussData:
    with open(file_path) as f:
        data = json.load(f)

    nodes = {}
    for node_data in data["nodes"]:
        constraints = node_data.get("constraints", "")
        constrained_x = "x" in constraints
        constrained_y = "y" in constraints

        deformations = node_data.get("deformations", {})
        deformation_x = deformations.get("x", 0.0)
        deformation_y = deformations.get("y", 0.0)

        loads = node_data.get("loads", {})
        load_x = loads.get("x", 0.0)
        load_y = loads.get("y", 0.0)

        node = Node(
            id=node_data["id"],
            dx=eval(node_data["dx"]),
            dy=eval(node_data["dy"]),
            constrained_x=constrained_x,
            constrained_y=constrained_y,
            deformation_x=deformation_x,
            deformation_y=deformation_y,
            load_x=load_x,
            load_y=load_y,
        )

        nodes[node.id] = node


    elements = [
        Element(
            starting_node=element["starting_node"], ending_node=element["ending_node"]
        )
        for element in data["elements"]
    ]
    return TrussData(
        nodes,
        elements
    )


truss: TrussData = load_data("./data/input.json")

render_truss_structure(truss)

solver = TrussSolver(truss)

solver.solve()

""" print(
    f"First element has cos: {truss.elements[0].cos} and sin: {truss.elements[0].sin}"
)
print("Stiffness matrix of the first element:")
print(truss.elements[0].stiffness_matrix) """
