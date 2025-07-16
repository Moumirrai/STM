from dataclasses import dataclass
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
    dof_x: Optional[int] = None
    dof_y: Optional[int] = None
    local_deformations: Optional[np.ndarray] = None


@dataclass
class Element:
    nodes: tuple[Node, Node]
    E: float = 210e6
    A: float = 0.01
    local_deformations: Optional[np.ndarray] = None

    def DOFs(self) -> List[int]:
        dofs = []
        for node in self.nodes:
            node_dofs = [node.dof_x, node.dof_y]
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
        base = (self.E * self.A) / self.magnitude()

        cos, sin = self.get_cos_sin()

        uu = base * cos ** 2
        uw = base * cos * sin
        ww = base * sin ** 2

        local_matrix = np.array([[uu, uw], [uw, ww]])

        return np.block([[local_matrix, -local_matrix], [-local_matrix, local_matrix]])

    def set_local_deformations(self, deformations: csc_array) -> None:

        self.local_deformations = np.zeros(len(self.nodes) * 2, dtype=float)

        for i, node in enumerate(self.nodes):
            node_deformations = np.zeros(2, dtype=float)
            if node.dof_x < deformations.shape[0]:  # if dof is free
                self.local_deformations[i * 2] = deformations[node.dof_x]
                node_deformations[0] = deformations[node.dof_x]
            if node.dof_y < deformations.shape[0]:
                self.local_deformations[i * 2 + 1] = deformations[node.dof_y]
                node_deformations[1] = deformations[node.dof_y]
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
