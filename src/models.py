from dataclasses import dataclass, field
from typing import Dict, List, Optional
from scipy.sparse import lil_matrix, linalg
from scipy.sparse import csc_array
import numpy as np


@dataclass
class Node:
    index: int
    dx: float
    dy: float
    constrained_x: bool = False
    constrained_y: bool = False
    deformation_x: float = 0.0
    deformation_y: float = 0.0
    load_x: float = 0.0
    load_y: float = 0.0
    dependency: Optional['Dependency'] = None
    local_deformations: Optional[np.ndarray] = None
    eigenstrain: np.ndarray = field(default_factory=lambda: np.zeros(2, dtype=float))


@dataclass
class MasterNode:
    nodeIndex: int
    direction: int # 0 for x, 1 for y
    factor: float
    eigenstrain: bool = True

@dataclass
class Dependency:
    dependency_index: int
    dependant_x: bool
    dependant_y: bool
    masters: List[MasterNode]


@dataclass
class Element:
    nodes: tuple[Node, Node]
    E: float = 210e6
    A: float = 0.01
    local_deformations: Optional[np.ndarray] = None
    _stiffness_matrix: Optional[np.ndarray] = None

    def __post_init__(self):
        displacements = np.array([
            self.nodes[0].deformation_x,
            self.nodes[0].deformation_y,
            self.nodes[1].deformation_x,
            self.nodes[1].deformation_y
        ])

        if not np.all(displacements == 0):
            forces = self.stiffness() @ displacements

            self.nodes[0].load_x += forces[0]
            self.nodes[0].load_y += forces[1]
            self.nodes[1].load_x += forces[2]
            self.nodes[1].load_y += forces[3]

        print(f"Initial forces for nodes {self.nodes[0].index} and {self.nodes[1].index} are {self.nodes[0].load_x}, {self.nodes[0].load_y}, {self.nodes[1].load_x}, {self.nodes[1].load_y}")

    def getDOFs(self) -> List[int]:
        dofs = []
        for node in self.nodes:
            # since dof ids are in same order as nodes we can get both from node index
            node_dofs = [node.index * 2, node.index * 2 + 1]
            dofs.extend(node_dofs)
        return dofs

    def magnitude(self) -> float:
        return (
                (self.nodes[0].dx - self.nodes[1].dx) ** 2
                + (self.nodes[0].dy - self.nodes[1].dy) ** 2
        ) ** 0.5

    def get_cos_sin(self) -> tuple[float, float]:
        vector_ab = (
            self.nodes[1].dx - self.nodes[0].dx,
            self.nodes[1].dy - self.nodes[0].dy,
        )

        length = self.magnitude()
        cos = vector_ab[0] / length
        sin = vector_ab[1] / length

        return cos, sin

    def stiffness(
            self,
    ) -> np.ndarray:
        if self._stiffness_matrix is not None:
            return self._stiffness_matrix
        base = (self.E * self.A) / self.magnitude()

        cos, sin = self.get_cos_sin()

        uu = base * cos ** 2
        uw = base * cos * sin
        ww = base * sin ** 2

        local_matrix = np.array([[uu, uw], [uw, ww]])

        self._stiffness_matrix = np.block([[local_matrix, -local_matrix], [-local_matrix, local_matrix]])
        return self._stiffness_matrix

    def set_local_deformations(self, deformations: csc_array) -> None:

        self.local_deformations = np.zeros(len(self.nodes) * 2, dtype=float)

        for i, node in enumerate(self.nodes):
            node_deformations = np.zeros(2, dtype=float)
            dof_x = node.index * 2
            dof_y = dof_x + 1
            if dof_x < deformations.shape[0]:  # if dof is free
                self.local_deformations[i * 2] = deformations[dof_x]
                node_deformations[0] = deformations[dof_x]
            if dof_y < deformations.shape[0]:
                self.local_deformations[i * 2 + 1] = deformations[dof_y]
                node_deformations[1] = deformations[dof_y]
            node.local_deformations = node_deformations

    def axial_force(self) -> float:
        cos, sin = self.get_cos_sin()
        stiffness_matrix = self.stiffness()
        forces = stiffness_matrix @ self.local_deformations
        normal_force = np.dot(forces[2:], np.array([cos, sin]))
        return normal_force


@dataclass
class TrussData:
    nodes: List[Node]
    elements: List[Element]
    constrained_dofs_count: int
    eigenstrainVector: Optional[np.ndarray] = (0, 0, 0)
