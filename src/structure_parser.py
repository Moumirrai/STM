import json
from typing import List, Optional

import numpy as np
import math

from models import TrussData, Node, Dependency, MasterNode, Element 


def parse_json_file(file_path: str, explicitEigenStrain: Optional[np.ndarray] = None) -> TrussData:
    with open(file_path) as f:
        data = json.load(f)

    total_constraints = 0
    
    default_E = data.get("defaultYoungsModulus", 210e6)
    default_A = data.get("defaultCrossSectionArea", 0.000004)

    nodes: List[Node] = []
    for i, node_data in enumerate(data["nodes"]):
        constraints = node_data.get("constraints", "")
        constrained_x = "x" in constraints
        constrained_y = "y" in constraints

        total_constraints += int(constrained_x) + int(constrained_y)

        deformations = node_data.get("deformations", {})
        deformation_x = deformations.get("x", 0.0)
        deformation_y = deformations.get("y", 0.0)

        loads = node_data.get("loads", {})
        load_x = loads.get("x", 0.0)
        load_y = loads.get("y", 0.0)

        new_node = Node(
            index=i,
            dx=eval(node_data["dx"]), # using eval to allow fractions
            dy=eval(node_data["dy"]),
            constrained_x=constrained_x,
            constrained_y=constrained_y,
            deformation_x=deformation_x,
            deformation_y=deformation_y,
            load_x=load_x,
            load_y=load_y,
            eigenstrain=np.array([0.0, 0.0]),
        )

        nodes.append(new_node)

    if explicitEigenStrain is not None:
        eigenstrain_vector = explicitEigenStrain
    else:
        eigenstrain_data = data.get("eigenstrain", {})
        eigenstrain_vector = np.array([
            eigenstrain_data.get("x", 0.0),
            eigenstrain_data.get("y", 0.0),
            eigenstrain_data.get("angle", 0.0)
        ])

    for i, dependency in enumerate(data.get("dependencies", [])):
        node_index = dependency["node"]
        node = nodes[node_index]
        
        if len(dependency["masters"]) == 0:
            print(f"Warning: Dependency {i} for node {node_index} has no masters.")
            continue

        if node.dependency is None:
            node.dependency = Dependency(
                masters=[],
                dependant_x=False,
                dependant_y=False,
                dependency_index=i
            )

        master_nodes = []
        for master_data in dependency["masters"]:
            direction = 0 if master_data["direction"] == "x" else 1

            if direction == 0:
                node.dependency.dependant_x = True
            else:
                node.dependency.dependant_y = True

            if master_data.get("eigenstrain") != False: # this is correct, we default to True if not specified
                master_node = nodes[master_data["node"]]

                if direction == 0:  # x direction
                    node.eigenstrain[0] += -1 * ((master_node.dx - node.dx) * eigenstrain_vector[0] + math.tan(eigenstrain_vector[2])/2 * (master_node.dy - node.dy))
                else:  # y direction
                    node.eigenstrain[1] += -1 * ((master_node.dy - node.dy) * eigenstrain_vector[1] + math.tan(eigenstrain_vector[2])/2 * (master_node.dx - node.dx))

            master_nodes.append(MasterNode(
                nodeIndex=master_data["node"],
                factor=master_data["factor"],
                direction=direction
            ))

        node.dependency.masters.extend(master_nodes)

    #volume = data.get("volume", 1.0) # TODO: calculate volume from gemotry if not provided
    max_x = max(node.dx for node in nodes)
    max_y = max(node.dy for node in nodes)
    min_x = min(node.dx for node in nodes)
    min_y = min(node.dy for node in nodes)
    volume = (max_x - min_x) * (max_y - min_y)

    elements = [
        Element(
            nodes=(
                nodes[int(element["starting_node"])],
                nodes[int(element["ending_node"])],
            ),
            E=element.get("E", default_E),
            A=element.get("A", default_A)
        )
        for element in data["elements"]
    ]

    return TrussData(nodes, elements, total_constraints, volume)