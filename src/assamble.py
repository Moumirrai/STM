# read data from nodes.inp
import math
from pathlib import Path

import numpy as np

from models import TrussData
from structure_parser import DependencyDefinition, EigenstrainDefinition, ElementDefinition, MasterDefinition, NodeDefinition, StructureDefinition, parse_structure_data, truss_data_to_json_string

def read_coords(file_path):
    coords = []
    with open(file_path, 'r') as f:
        lines = f.readlines()[1:]
        for line in lines:
            parts = line.strip().split('\t')[1:3]
            coords.append([float(part) for part in parts])
    return coords

def parse_elements(file_path):
    elements = []
    with open(file_path, 'r') as f:
        lines = f.readlines()[1:]
        for line in lines:
            parts = line.strip().split('\t')[1:3] + line.strip().split('\t')[4:6]
            elements.append([int(part) for part in parts])
    return elements

def parse_dependencies(file_path):
    dependencies = []
    with open(file_path, 'r') as f:
        line = f.readline()
        parts = line.strip().split('\t')[14:]
        values = [int(part) for part in parts]
        if len(values) % 2 != 0:
            raise ValueError("Expected an even number of dependency indices")
        dependencies.extend([values[i:i + 2] for i in range(0, len(values), 2)])
    return dependencies

raw_nodes = read_coords('ref_data/nodes.inp')
auxNodes = read_coords('ref_data/auxNodes.inp')
vertices = read_coords('ref_data/vertices.inp')
raw_dependencies = parse_dependencies('ref_data/blocks.inp')

#get width and height of the structure from the vertices
width = max(vertex[0] for vertex in vertices) - min(vertex[0] for vertex in vertices)
height = max(vertex[1] for vertex in vertices) - min(vertex[1] for vertex in vertices)

coordList = raw_nodes + auxNodes + vertices

raw_elements = parse_elements('ref_data/mechElems.inp')

elements = []
for element in raw_elements:
    x1, y1 = coordList[element[2]]
    x2, y2 = coordList[element[3]]
    length = math.hypot(x2 - x1, y2 - y1)
    elements.append(element[:2] + [length])
    
node_definitions = [NodeDefinition(dx=x, dy=y) for x, y in raw_nodes]
node_definitions[0].constraints = "xy"

element_definitions = []
for index_a, index_b, length in elements:
    # Cross-section proportional to Voronoi ridge length; E is set globally via defaultYoungsModulus
    element_definitions.append(ElementDefinition(
        starting_node=index_a,
        ending_node=index_b,
        A=length,
    ))
    
dependencies = [DependencyDefinition(
    node=ghost_index,
    masters=[MasterDefinition(node=master_index, direction="x", factor=1.0), MasterDefinition(node=master_index, direction="y", factor=1.0)]
) for master_index, ghost_index in raw_dependencies]


    
structure = StructureDefinition(
    nodes=node_definitions,
    elements=element_definitions,
    dependencies=dependencies,
    eigenstrain=EigenstrainDefinition(x=1.0, y=0.0, angle=1.0),
    defaultYoungsModulus=210e9,
    volume=width * height,
)

truss: TrussData = parse_structure_data(
    structure, explicitEigenStrain=np.array([0.0, 0.0, 0.0])
)

jsonString = truss_data_to_json_string(truss)

output_path = Path(__file__).resolve().parent.parent / "structures" / "-1.json"
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text(jsonString, encoding="utf-8")