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
        """Solve the truss system"""

        total_dof_count = len(self.truss.nodes) * 2

        self.DOFIDS = np.zeros(total_dof_count, dtype=int)
        free_dof_counter = 0
        total_free_dof_count = total_dof_count - self.truss.constrained_dofs_count

        fixed_dof_counter = total_free_dof_count # no need to add 1 since free dof counter starts from 0

        S = np.zeros((total_free_dof_count), dtype=float)

        r = np.zeros((total_free_dof_count), dtype=float)

        # create a mapping of dof IDs
        for i, node in enumerate(self.truss.nodes):
            node_dof_id = i * 2 

            if not node.constrained_x:
                node.dof_x = free_dof_counter
                self.DOFIDS[node_dof_id] = free_dof_counter
                if node.load_x:
                    S[free_dof_counter] = node.load_x
                if node.deformation_x:
                    r[free_dof_counter] = node.deformation_x
                free_dof_counter += 1
            else:
                node.dof_x = fixed_dof_counter

                self.DOFIDS[node_dof_id] = fixed_dof_counter
                fixed_dof_counter += 1
                

            if not node.constrained_y:
                node.dof_y = free_dof_counter
                self.DOFIDS[node_dof_id + 1] = (
                    free_dof_counter
                )
                if node.load_y:
                    S[free_dof_counter] = node.load_y
                if node.deformation_y:
                    r[free_dof_counter] = node.deformation_y
                free_dof_counter += 1
            else:
                node.dof_y = fixed_dof_counter
                self.DOFIDS[node_dof_id + 1] = (
                    fixed_dof_counter
                )
                fixed_dof_counter += 1

            

        # we well use the row-based list of lists sparse matrix for initial construction
        # and then convert it to csr format for more efficient operations
        lil_K = lil_matrix((free_dof_counter, free_dof_counter), dtype=float)

        for element in self.truss.elements:
            stiffness_matrix = element.stiffness()

            dofs = element.DOFs()

            # populate the global stiffness matrix
            for i, row_dof in enumerate(dofs):
                if row_dof < free_dof_counter:  # skip constrained dofs
                    for j, col_dof in enumerate(dofs):
                        if col_dof < free_dof_counter:  # skip constrained dofs
                            lil_K[row_dof, col_dof] += stiffness_matrix[i, j]

        #print("S", S)

        self.K = (
            lil_K.tocsr()
        )  # convert to compressed sparse row format for more efficient operations

        u_free = linalg.spsolve(self.K, S)

        for element in self.truss.elements:
            element.set_local_deformations(u_free)

        #self.check_stability()

        #print(self.K.toarray())

        pass
