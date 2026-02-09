from typing import List, Dict

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from dataclasses import dataclass
from utils import dump_matrix_to_csv

from models import TrussData

@dataclass
class FixedConstraint:
    dof: int
    value: float = 0.0

@dataclass
class DependentConstraint:
    dof: int
    coeffs: Dict[int, float]
    eigenstrain: np.ndarray

class LagrangeTrussSolver:
    def __init__(self, truss: TrussData):
        self.truss = truss

    def solve(self) -> np.ndarray:

        total_dof_count = len(self.truss.nodes) * 2
        fixedConstraints: List[FixedConstraint] = []
        dependentConstraints: List[DependentConstraint] = []
        f = np.zeros(total_dof_count)

        for node in self.truss.nodes:
            base_dof = node.index * 2
            if node.constrained_x:
                fixedConstraints.append(FixedConstraint(dof=base_dof, value=node.deformation_x))
            if node.constrained_y:
                fixedConstraints.append(FixedConstraint(dof=base_dof + 1, value=node.deformation_y))
            # f vector
            f[base_dof] = node.load_x
            f[base_dof + 1] = node.load_y

        #print(f)

        for node in self.truss.nodes:
            base_dof = node.index * 2
            if node.dependency:
                # X direction
                if node.dependency.dependant_x:
                    x_dof = base_dof
                    coeffs = {}
                    eigenstrain = node.eigenstrain[0]  # From node.eigenstrain (assumed [x, y])
                    for master in node.dependency.masters:
                        if master.direction == 0:  # X direction
                            master_dof = master.nodeIndex * 2
                            coeffs[master_dof] = master.factor
                    dependentConstraints.append(DependentConstraint(
                        dof=x_dof,
                        eigenstrain=eigenstrain,
                        coeffs=coeffs,
                    ))

                # Y direction
                if node.dependency.dependant_y:
                    y_dof = base_dof + 1
                    coeffs = {}
                    eigenstrain = node.eigenstrain[1]
                    for master in node.dependency.masters:
                        if master.direction == 1:  # Y direction
                            master_dof = master.nodeIndex * 2 + 1
                            coeffs[master_dof] = master.factor
                    dependentConstraints.append(DependentConstraint(
                        dof=y_dof,
                        eigenstrain=eigenstrain,
                        coeffs=coeffs,
                    ))

        num_constraints = len(fixedConstraints) + len(dependentConstraints)

        g = np.zeros(num_constraints)

        K = lil_matrix((total_dof_count, total_dof_count))
        for element in self.truss.elements:
            stiffness_matrix = element.stiffness()  # 4x4 matrix from element
            dofs = element.getDOFs()  # List of 4 global DOF indices
            for i in range(4):
                for j in range(4):
                    K[dofs[i], dofs[j]] += stiffness_matrix[i, j]  # Accumulate (handles overlapping elements)
        K = K.tocsr()  # Convert to CSR for efficient operations

        C = lil_matrix((num_constraints, total_dof_count))
        row = 0
        # Fixed constraints
        for constraint in fixedConstraints:
            C[row, constraint.dof] = 1
            g[row] = constraint.value
            row += 1
        # Dependent constraints (keep existing loop, but update row)
        for constraint in dependentConstraints:
            C[row, constraint.dof] = 1
            for master_dof, factor in constraint.coeffs.items():
                C[row, master_dof] = -factor  # Fix: negative for masters
            g[row] = constraint.eigenstrain
            row += 1

        C = C.tocsr()
        f_aug = np.concatenate([f, g])

        aug_size = total_dof_count + num_constraints
        K_aug = lil_matrix((aug_size, aug_size))
        K_aug[:total_dof_count, :total_dof_count] = K
        K_aug[:total_dof_count, total_dof_count:] = C.T  # Remove * -1
        K_aug[total_dof_count:, :total_dof_count] = C
        K_aug = K_aug.tocsr()

        u_aug = spsolve(K_aug, f_aug)

        u_vec_solved = u_aug[:total_dof_count]

        lambdas = u_aug[total_dof_count:]

        #dump_matrix_to_csv(K_aug, "debug_export.csv")

        stress_contributions = []

        for element in self.truss.elements:
            element.set_local_deformations(u_vec_solved)

            forces = element.forces_vector()

            value = (
                element.magnitude()
                * element.axial_force()
                * np.multiply.outer(element.get_cos_sin(), element.get_cos_sin())
            )
            stress_contributions.append(value)

        result = 1 / self.truss.volume * sum(stress_contributions)
        # print("Result:")
        # print(result)

        # Extract stress components (xx, yy, xy) from the 2x2 result matrix
        return np.array(
            [
                result[0, 0],  # xx
                result[1, 1],  # yy
                result[0, 1] * 2,  # xy
            ]
        )
