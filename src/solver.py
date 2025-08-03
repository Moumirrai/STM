import numpy as np
from models import Element, TrussData
from scipy.sparse import lil_matrix, linalg, identity, coo_matrix

E = 210e6  # Pa
A = 0.01  # m^2


class TrussSolver:

    def __init__(
            self, truss: TrussData, young_modulus: float = 210e6, area: float = 0.01
    ):
        self.truss = truss
        self.K = None  # global stiffness matrix
        self.DOFIDS = None  # degrees of freedom IDs

    def check_stability(self):
        """Check if the truss is stable"""
        if self.K is None:
            raise ValueError("Stiffness matrix is not computed yet.")

        eigenvalues, _ = linalg.eigs(
            self.K, k=1, which="SM"
        )  # check for smallest eigenvalue
        min_eigenvalue = eigenvalues.real[0]
        if min_eigenvalue <= 1e-9:
            raise ValueError(
                "Truss is unstable! The smallest eigenvalue is close to zero."
            )

    def solve(self):

        total_dof_count = len(self.truss.nodes) * 2

        dependency_map = {dep.nodeIndex: dep for dep in self.truss.dependencies}

        # arrays to hold indices of free, dependent, and fixed dofs
        free_dof_indices = []
        dependent_dof_indices = []
        fixed_dof_indices = []

        for node_idx, node in enumerate(self.truss.nodes):
            base_dof = node_idx * 2

            # X direction DOF
            x_dof = base_dof
            if node_idx in dependency_map:
                dependent_dof_indices.append(x_dof)
            elif node.constrained_x:
                fixed_dof_indices.append(x_dof)
            else:
                free_dof_indices.append(x_dof)

            # Y direction DOF
            y_dof = base_dof + 1
            if node_idx in dependency_map:
                dependent_dof_indices.append(y_dof)
            elif node.constrained_y:
                fixed_dof_indices.append(y_dof)
            else:
                free_dof_indices.append(y_dof)

        # concat all indices
        ordered_dof_indices = free_dof_indices + dependent_dof_indices + fixed_dof_indices

        print(f"Free DOFs: {len(free_dof_indices)}, Dependent DOFs: {len(dependent_dof_indices)}, Fixed DOFs: {len(fixed_dof_indices)}")
        print(f"Ordered DOF indices: {ordered_dof_indices}")
            

