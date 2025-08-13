import json
from typing import List

from models import TrussData, Node, Dependency, MasterNode, Element


def parse_json_file(file_path: str) -> TrussData:
    with open(file_path) as f:
        data = json.load(f)

    total_constraints = 0

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
            dx=eval(node_data["dx"]),
            dy=eval(node_data["dy"]),
            constrained_x=constrained_x,
            constrained_y=constrained_y,
            deformation_x=deformation_x,
            deformation_y=deformation_y,
            load_x=load_x,
            load_y=load_y,
        )

        nodes.append(new_node)

    for i, dependency in enumerate(data.get("dependencies", [])):
        node_index = dependency["node"]
        node = nodes[node_index]

        dependency_object = Dependency(
            masters= [],
            dependant_x=False,
            dependant_y=False,
            dependency_index=i
        )

        if len(dependency["masters"]) == 0:
            print(f"Warning: Dependency {i} for node {node_index} has no masters.")
            continue

        # Process masters and set dependency flags in one pass
        master_nodes = []
        for master_data in dependency["masters"]:
            direction = 0 if master_data["direction"] == "x" else 1

            # Set dependency flags based on direction
            if direction == 0:
                dependency_object.dependant_x = True
            else:
                dependency_object.dependant_y = True

            master_nodes.append(MasterNode(
                nodeIndex=master_data["node"],
                factor=master_data["factor"],
                direction=direction
            ))

        dependency_object.masters = master_nodes

        node.dependency = dependency_object

    elements = [
        Element(
            nodes=(
                nodes[int(element["starting_node"])],
                nodes[int(element["ending_node"])],
            ),
        )
        for element in data["elements"]
    ]
    return TrussData(nodes, elements, total_constraints)