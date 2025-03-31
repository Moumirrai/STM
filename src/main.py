import json
from dataclasses import dataclass
from typing import Dict, List, Optional
import matplotlib.pyplot as plt
import numpy as np

@dataclass
class Node:
    id: int
    dx: float
    dy: float

@dataclass
class Element:
    starting_node: int
    ending_node: int
    sin: Optional[float] = None
    cos: Optional[float] = None
    stiffness_matrix: Optional[np.ndarray] = None
    length: Optional[float] = None

@dataclass
class TrussData:
    nodes: Dict[int, Node]
    elements: List[Element]

def load_data(file_path: str) -> TrussData:
    with open(file_path) as f:
        data = json.load(f)

    nodes = {
        int(node['id']): Node(
            id=int(node['id']),
            dx=eval(node['dx']),
            dy=eval(node['dy'])
        )
        for node in data['nodes']
    }
    elements = [
        Element(
            starting_node=element['starting_node'],
            ending_node=element['ending_node']
        )
        for element in data['elements']
    ]
    return TrussData(nodes, elements)

truss: TrussData = load_data("./data/input.json")


### PLOT
for node in truss.nodes.values():
    plt.plot(node.dx, node.dy, 'o', color='black')
    plt.text(node.dx, node.dy - 0.02, f"N{node.id}", fontsize=8, ha='right')

for element in truss.elements:
    start_node = truss.nodes[element.starting_node]
    end_node = truss.nodes[element.ending_node]
    plt.plot([start_node.dx, end_node.dx], [start_node.dy, end_node.dy], 'k-')

plt.xlabel('dx')
plt.ylabel('dy')
plt.title('Nodes and Elements')
plt.grid(True)

plt.gca().invert_yaxis()

#plt.show()


### END PLOT

def compute_angle(element: Element):
    start_node = truss.nodes[element.starting_node]
    end_node = truss.nodes[element.ending_node]

    vectorAB = (end_node.dx - start_node.dx, end_node.dy - start_node.dy)

    magnitudeAB = (vectorAB[0]**2 + vectorAB[1]**2)**0.5

    element.length = magnitudeAB

    element.cos = vectorAB[0] / magnitudeAB
    element.sin = vectorAB[1] / magnitudeAB

E = 210e6  # Pa
A = 0.01  # m^2

def compute_stiffness(element: Element):

    base = (E * A)/element.length

    uu = base * element.cos ** 2
    uw = base * element.cos * element.sin
    ww = base * element.sin ** 2


    local_matrix = np.array([[uu, uw],
                             [uw, ww]])
    
    element.stiffness_matrix = np.block([
        [ local_matrix, -local_matrix ],
        [ -local_matrix, local_matrix  ]
    ])

for element in truss.elements:
    compute_angle(element)
    compute_stiffness(element)


print(f"First element has cos: {truss.elements[0].cos} and sin: {truss.elements[0].sin}")
print("Stiffness matrix of the first element:")
print(truss.elements[0].stiffness_matrix)

plt.show()