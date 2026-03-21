import json
from typing import List, Optional, Tuple, Dict, Any
from dataclasses import dataclass, field

import numpy as np
import math

from models import TrussData, Node, Dependency, MasterNode, Element 


@dataclass
class NodeDefinition:
    dx: float
    dy: float
    constraints: str = ""
    deformations: Dict[str, float] = field(default_factory=dict)
    loads: Dict[str, float] = field(default_factory=dict)


@dataclass
class ElementDefinition:
    starting_node: int
    ending_node: int
    E: Optional[float] = None
    A: Optional[float] = None


@dataclass
class MasterDefinition:
    node: int
    direction: str  # "x" or "y"
    factor: float
    eigenstrain: bool = True


@dataclass
class DependencyDefinition:
    node: int
    masters: List[MasterDefinition]


@dataclass
class EigenstrainDefinition:
    x: float = 0.0
    y: float = 0.0
    angle: float = 0.0


@dataclass
class StructureDefinition:
    nodes: List[NodeDefinition]
    elements: List[ElementDefinition]
    dependencies: List[DependencyDefinition] = field(default_factory=list)
    eigenstrain: EigenstrainDefinition = field(default_factory=EigenstrainDefinition)
    defaultYoungsModulus: float = 210e6
    defaultCrossSectionArea: float = 0.000004
    volume: Optional[float] = None

    @classmethod
    def from_json_dict(cls, data: dict) -> 'StructureDefinition':
        """Create StructureDefinition from JSON-like dict format."""
        nodes = []
        for node_data in data.get("nodes", []):
            nodes.append(NodeDefinition(
                dx=eval(node_data["dx"]),  # Still need eval for fractions
                dy=eval(node_data["dy"]),
                constraints=node_data.get("constraints", ""),
                deformations=node_data.get("deformations", {}),
                loads=node_data.get("loads", {})
            ))
        
        elements = []
        for elem_data in data.get("elements", []):
            elements.append(ElementDefinition(
                starting_node=int(elem_data["starting_node"]),
                ending_node=int(elem_data["ending_node"]),
                E=elem_data.get("E"),
                A=elem_data.get("A")
            ))
        
        dependencies = []
        for dep_data in data.get("dependencies", []):
            masters = []
            for master_data in dep_data.get("masters", []):
                masters.append(MasterDefinition(
                    node=master_data["node"],
                    direction=master_data["direction"],
                    factor=master_data["factor"],
                    eigenstrain=master_data.get("eigenstrain", True)
                ))
            dependencies.append(DependencyDefinition(
                node=dep_data["node"],
                masters=masters
            ))
        
        eigenstrain_data = data.get("eigenstrain", {})
        eigenstrain = EigenstrainDefinition(
            x=eigenstrain_data.get("x", 0.0),
            y=eigenstrain_data.get("y", 0.0),
            angle=eigenstrain_data.get("angle", 0.0)
        )
        
        return cls(
            nodes=nodes,
            elements=elements,
            dependencies=dependencies,
            eigenstrain=eigenstrain,
            defaultYoungsModulus=data.get("defaultYoungsModulus", 210e6),
            defaultCrossSectionArea=data.get("defaultCrossSectionArea", 0.000004),
            volume=data.get("volume")
        )

    def to_truss_data(self, explicitEigenStrain: Optional[np.ndarray] = None) -> TrussData:
        """Parse this structure definition into a TrussData object."""
        return parse_structure_data(self, explicitEigenStrain) 


def parse_structure_data(definition: StructureDefinition, explicitEigenStrain: Optional[np.ndarray] = None) -> TrussData:
    total_constraints = 0
    
    default_E = definition.defaultYoungsModulus
    default_A = definition.defaultCrossSectionArea

    nodes: List[Node] = []
    for i, node_def in enumerate(definition.nodes):
        constraints = node_def.constraints
        constrained_x = "x" in constraints
        constrained_y = "y" in constraints

        total_constraints += int(constrained_x) + int(constrained_y)

        deformations = node_def.deformations
        deformation_x = deformations.get("x", 0.0)
        deformation_y = deformations.get("y", 0.0)

        loads = node_def.loads
        load_x = loads.get("x", 0.0)
        load_y = loads.get("y", 0.0)

        new_node = Node(
            index=i,
            dx=node_def.dx,
            dy=node_def.dy,
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
        eigenstrain_vector = np.array([
            definition.eigenstrain.x,
            definition.eigenstrain.y,
            definition.eigenstrain.angle
        ])

    for dep_def in definition.dependencies:
        node_index = dep_def.node
        node = nodes[node_index]
        
        if len(dep_def.masters) == 0:
            print(f"Warning: Dependency for node {node_index} has no masters.")
            continue

        if node.dependency is None:
            node.dependency = Dependency(
                masters=[],
                dependant_x=False,
                dependant_y=False,
                dependency_index=len(definition.dependencies)  # Use index in list
            )

        for master_def in dep_def.masters:
            direction = 0 if master_def.direction == "x" else 1

            if direction == 0:
                node.dependency.dependant_x = True
            else:
                node.dependency.dependant_y = True

            if master_def.eigenstrain:
                master_node = nodes[master_def.node]

                if direction == 0:  # x direction
                    node.eigenstrain[0] += -1 * ((master_node.dx - node.dx) * eigenstrain_vector[0] + math.tan(eigenstrain_vector[2])/2 * (master_node.dy - node.dy))
                else:  # y direction
                    node.eigenstrain[1] += -1 * ((master_node.dy - node.dy) * eigenstrain_vector[1] + math.tan(eigenstrain_vector[2])/2 * (master_node.dx - node.dx))

            master_nodes = MasterNode(
                nodeIndex=master_def.node,
                factor=master_def.factor,
                direction=direction
            )

            node.dependency.masters.append(master_nodes)

    # Calculate volume if not provided
    if definition.volume is not None:
        volume = definition.volume
    else:
        max_x = max(node.dx for node in nodes)
        max_y = max(node.dy for node in nodes)
        min_x = min(node.dx for node in nodes)
        min_y = min(node.dy for node in nodes)
        volume = (max_x - min_x) * (max_y - min_y)

    elements = [
        Element(
            nodes=(
                nodes[elem_def.starting_node],
                nodes[elem_def.ending_node],
            ),
            E=elem_def.E if elem_def.E is not None else default_E,
            A=elem_def.A if elem_def.A is not None else default_A
        )
        for elem_def in definition.elements
    ]

    return TrussData(nodes, elements, total_constraints, volume)


def parse_json_file(file_path: str, explicitEigenStrain: Optional[np.ndarray] = None) -> TrussData:
    with open(file_path) as f:
        data = json.load(f)
    definition = StructureDefinition.from_json_dict(data)
    return parse_structure_data(definition, explicitEigenStrain)


def truss_data_to_json_string(truss: TrussData, indent: Optional[int] = 2) -> str:
    """Serialize TrussData into a JSON string compatible with StructureDefinition-style input."""

    def _is_close(a: float, b: float) -> bool:
        return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-15)

    def _format_constraints(node: Node) -> str:
        constraints = ""
        if node.constrained_x:
            constraints += "x"
        if node.constrained_y:
            constraints += "y"
        return constraints

    # Infer defaults from elements when possible so JSON stays concise.
    default_E: Optional[float] = None
    default_A: Optional[float] = None
    if truss.elements:
        first_E = truss.elements[0].E
        first_A = truss.elements[0].A
        if all(_is_close(elem.E, first_E) for elem in truss.elements):
            default_E = first_E
        if all(_is_close(elem.A, first_A) for elem in truss.elements):
            default_A = first_A

    nodes_data: List[Dict[str, Any]] = []
    for node in truss.nodes:
        node_data: Dict[str, Any] = {
            "dx": str(node.dx),
            "dy": str(node.dy),
        }

        constraints = _format_constraints(node)
        if constraints:
            node_data["constraints"] = constraints

        deformations: Dict[str, float] = {}
        if node.deformation_x != 0.0:
            deformations["x"] = node.deformation_x
        if node.deformation_y != 0.0:
            deformations["y"] = node.deformation_y
        if deformations:
            node_data["deformations"] = deformations

        loads: Dict[str, float] = {}
        if node.load_x != 0.0:
            loads["x"] = node.load_x
        if node.load_y != 0.0:
            loads["y"] = node.load_y
        if loads:
            node_data["loads"] = loads

        nodes_data.append(node_data)

    elements_data: List[Dict[str, Any]] = []
    for element in truss.elements:
        element_data: Dict[str, Any] = {
            "starting_node": element.nodes[0].index,
            "ending_node": element.nodes[1].index,
        }

        if default_E is None or not _is_close(element.E, default_E):
            element_data["E"] = element.E
        if default_A is None or not _is_close(element.A, default_A):
            element_data["A"] = element.A

        elements_data.append(element_data)

    dependencies_data: List[Dict[str, Any]] = []
    for node in truss.nodes:
        if node.dependency is None:
            continue

        masters_data = []
        for master in node.dependency.masters:
            masters_data.append(
                {
                    "node": master.nodeIndex,
                    "direction": "x" if master.direction == 0 else "y",
                    "factor": master.factor,
                    "eigenstrain": master.eigenstrain,
                }
            )

        if masters_data:
            dependencies_data.append(
                {
                    "node": node.index,
                    "masters": masters_data,
                }
            )

    output: Dict[str, Any] = {
        "nodes": nodes_data,
        "elements": elements_data,
        "dependencies": dependencies_data,
        "volume": truss.volume,
    }

    if default_E is not None:
        output["defaultYoungsModulus"] = default_E
    if default_A is not None:
        output["defaultCrossSectionArea"] = default_A

    return json.dumps(output, indent=indent)