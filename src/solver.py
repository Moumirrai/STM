import numpy as np
from models import Element, TrussData
from scipy.sparse import lil_matrix, linalg

E = 210e6  # Pa
A = 0.01  # m^2


class TrussSolver:

    def __init__(
        self, truss: TrussData, young_modulus: float = 210e6, area: float = 0.01
    ):
        self.truss = truss
        self.E = young_modulus  # Pa
        self.A = area  # m2
        self.K = None  # global stiffness matrix
        self.DOFIDS = None  # degrees of freedom IDs

    def compute_angle(self, element: Element):
        start_node = self.truss.nodes[element.starting_node]
        end_node = self.truss.nodes[element.ending_node]

        vectorAB = (end_node.dx - start_node.dx, end_node.dy - start_node.dy)
        magnitudeAB = (vectorAB[0] ** 2 + vectorAB[1] ** 2) ** 0.5

        element.length = magnitudeAB
        element.cos = vectorAB[0] / magnitudeAB
        element.sin = vectorAB[1] / magnitudeAB

    def compute_stiffness(
        self, element: Element
    ) -> np.ndarray:  # could be moved as a static method of the Element class
        base = (self.E * self.A) / element.length

        uu = base * element.cos**2
        uw = base * element.cos * element.sin
        ww = base * element.sin**2

        local_matrix = np.array([[uu, uw], [uw, ww]])

        return np.block([[local_matrix, -local_matrix], [-local_matrix, local_matrix]])

    def prepare_elements(self):
        """Compute angles"""
        for element in self.truss.elements:
            self.compute_angle(element)

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
        """Solve the truss system"""
        self.prepare_elements()

        total_dof_count = len(self.truss.nodes) * 2

        self.DOFIDS = np.zeros(total_dof_count, dtype=int)
        free_dof_count = 0

        # create a mapping og DOF IDs
        for i, node in self.truss.nodes.items():
            node_dof_id = (
                i - 1
            ) * 2  # index of the U direction DOF index for given node

            if not node.constrained_x:
                free_dof_count += 1
                self.DOFIDS[node_dof_id] = free_dof_count

            if not node.constrained_y:
                free_dof_count += 1
                self.DOFIDS[node_dof_id + 1] = (
                    free_dof_count  # increment by one to get the V direction DOF index
                )

        if (total_dof_count - free_dof_count) < 3:
            raise ValueError(
                "Structure is unstable! At least 3 DOFs must be constrained."
            )
        if len(self.truss.elements) <= (
            2 * len(self.truss.nodes) - (total_dof_count - free_dof_count)
        ):
            raise ValueError("Structure is likely unstable!")

        # we well use the row-based list of Lists sparse matrix for initial construction
        # and then convert it to Compressed Sparse Row format for more efficient operations
        lil_K = lil_matrix((free_dof_count, free_dof_count), dtype=float)

        for element in self.truss.elements:
            stiffness_matrix = self.compute_stiffness(element)

            start_node_dof_x = (element.starting_node - 1) * 2
            start_node_dof_y = start_node_dof_x + 1
            end_node_dof_x = (element.ending_node - 1) * 2
            end_node_dof_y = end_node_dof_x + 1

            dofs = [
                self.DOFIDS[start_node_dof_x],  # start node x
                self.DOFIDS[start_node_dof_y],  # start node y
                self.DOFIDS[end_node_dof_x],  # end node x
                self.DOFIDS[end_node_dof_y],  # end node y
            ]

            # populate the global stiffness matrix
            for i in range(4):
                row_dof = dofs[i]
                if row_dof > 0:  # skip constrained DOFs with zero ID
                    for j in range(4):
                        col_dof = dofs[j]
                        if col_dof > 0:  # skip constrained DOFs with zero ID
                            lil_K[row_dof - 1, col_dof - 1] += stiffness_matrix[i, j]

        self.K = (
            lil_K.tocsr()
        )  # convert to compressed sparse row format for more efficient operations

        self.check_stability()

        print(self.DOFIDS)

        pass
