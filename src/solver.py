import numpy as np
from models import TrussData
from scipy.sparse import lil_matrix, linalg, identity, bmat
from numpy import set_printoptions, save

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

        set_printoptions(
            linewidth=250,
        )
        total_dof_count = len(self.truss.nodes) * 2

        # arrays to hold indices of free, dependent, and fixed dofs
        free_dof_indices = []
        dependent_dof_indices = []
        fixed_dof_indices = []

        # iterate over nodes and classify their DOFs
        """
        Classification example with 3 nodes (6 DOFs):
          Free DOFs:      [0, 3]         (Node 0-x, Node 1-y)
          Dependent DOFs: [1, 4]         (Node 0-y, Node 2-x)
          Fixed DOFs:     [2, 5]         (Node 1-x, Node 2-y)
        
          Global DOF indices: [0, 1, 2, 3, 4, 5]
          Reordered DOFIDs:   [0, 3, 1, 4, 2, 5]
                               ^  ^  ^  ^  ^  ^
                               Free  Dependent  Fixed
        """
        for node in self.truss.nodes:

            # X direction DOF
            x_dof = node.index * 2
            if node.dependency and node.dependency.dependant_x:
                dependent_dof_indices.append(x_dof)
            elif node.constrained_x:
                fixed_dof_indices.append(x_dof)
            else:
                free_dof_indices.append(x_dof)

            # Y direction DOF
            y_dof = x_dof + 1
            if node.dependency and node.dependency.dependant_y:
                dependent_dof_indices.append(y_dof)
            elif node.constrained_y:
                fixed_dof_indices.append(y_dof)
            else:
                free_dof_indices.append(y_dof)

        # concat all indices
        self.DOFIDS = free_dof_indices + dependent_dof_indices + fixed_dof_indices

        """
        X matrix structure:
          +--------+--------+
          │  x_11  │  None  │ <- Free DOFs
          +--------+--------+
          │  x_d1  │  x_d2  │ <- Dependent DOFs
          +--------+--------+
          │  None  │  x_22  │ <- Fixed DOFs
          +--------+--------+
            Free     Fixed
        """

        x_11 = identity(len(free_dof_indices))
        x_22 = identity(len(fixed_dof_indices))

        free_and_fixed_indices = free_dof_indices + fixed_dof_indices

        # create the xd1+xd2 matrix
        x_d = lil_matrix((len(dependent_dof_indices), len(free_and_fixed_indices)))
        # create index map beforehand, so we don't have to index array every time
        reduced_index_map = {dof: idx for idx, dof in enumerate(free_and_fixed_indices)}

        for local_i, global_i in enumerate(dependent_dof_indices):
            # we get the node ID and direction from the global index by undoing the mapping
            node_id = global_i // 2
            direction = global_i % 2

            node = self.truss.nodes[node_id]

            filtered_masters = [  # filter only the direction we are currently processing
                master for master in node.dependency.masters
                if master.direction == direction
            ]

            for master in filtered_masters:
                master_global_dof = master.nodeIndex * 2 + master.direction  # convert to global index
                if master_global_dof in reduced_index_map:
                    master_local_index = reduced_index_map[master_global_dof]
                    # fox y coordinates we can increment loop index since we already have ordered the DOFs
                    # for x coordinates we use the master_local_index (note that this is the local index in the reduced matrix)
                    x_d[local_i, master_local_index] = master.factor

        # convert to CSR format
        x_d = x_d.tocsr()
        # now we can split the x_d matrix into x_d1 (free DOFs) and x_d2 (fixed DOFs)
        x_d1 = x_d[:, :len(free_dof_indices)]
        x_d2 = x_d[:, len(free_dof_indices):]

        # and assemble the full X matrix
        x_mat = bmat([
            [x_11, None],
            [x_d1, x_d2],
            [None, x_22]
        ], format='csr')

        r_reduced = np.zeros(len(free_and_fixed_indices))
        f_vec = np.zeros(len(free_dof_indices) + len(dependent_dof_indices))

        free_index_map = {dof: idx for idx, dof in enumerate(free_dof_indices)}
        fixed_index_map = {dof: idx for idx, dof in enumerate(fixed_dof_indices)}
        dependent_index_map = {dof: idx for idx, dof in enumerate(dependent_dof_indices)}

        for node_idx, node in enumerate(self.truss.nodes):
            base_dof = node_idx * 2
            loads = [node.load_x, node.load_y]
            deformations = [node.deformation_x, node.deformation_y]

            for direction, (load, deformation) in enumerate(zip(loads, deformations)):
                global_dof = base_dof + direction

                if global_dof in free_index_map:
                    r_reduced[free_index_map[global_dof]] = deformation
                    f_vec[free_index_map[global_dof]] = load
                elif global_dof in dependent_index_map:
                    f_vec[len(free_dof_indices) + dependent_index_map[global_dof]] = load
                elif global_dof in fixed_index_map:
                    r_reduced[len(free_dof_indices) + fixed_index_map[global_dof]] = deformation

        r_vec = x_mat.dot(r_reduced)

        # split r_vec into free and fixed parts
        r_fixed = r_vec[len(free_dof_indices) + len(dependent_dof_indices):]

        # now we need to assamble the global stiffness matrix K
        # we will devide it into 3 parts:
        # K11 - free DOFs, K1D - dependent DOFs, K12 - fixed DOFs, only the horizontal part is needed

        self.K = lil_matrix((total_dof_count, total_dof_count))

        """
        K matrix structure:
          +-------+-------+-------+
          │  K11  │  K1D  │  K12  │  <- Free
          +-------+-------+-------+
          │  KD1  │  KDD  │  KD2  │  <- Dep.
          +-------+-------+-------+
          │  K21  │  K2D  │  K22  │  <- Fixed
          +-------+-------+-------+
            Free     Dep.   Fixed
        """

        for element in self.truss.elements:
            stiffness_matrix = element.stiffness()

            # Global DOF indices for this element
            dofs = element.getDOFs()

            # Add to global stiffness matrix
            for i in range(4):
                for j in range(4):
                    self.K[dofs[i], dofs[j]] += stiffness_matrix[i, j]

        self.K = self.K.tocsr()

        # Partition K matrix according to DOF ordering
        K11 = self.K[np.ix_(free_dof_indices, free_dof_indices)]
        K1D = self.K[np.ix_(free_dof_indices, dependent_dof_indices)]
        K12 = self.K[np.ix_(free_dof_indices, fixed_dof_indices)]

        A = K11 + K1D @ x_d1
        B = K12 + K1D @ x_d2

        # Calculate right-hand side vector
        rhs = f_vec[:len(free_dof_indices)] - B @ r_fixed

        # Solve for free DOF displacements
        from scipy.sparse.linalg import spsolve
        r_free_solved = spsolve(A, rhs)

        print(f"Solved free displacements: {r_free_solved}")

        # Update the full displacement vector
        r_vec_solved = np.zeros(total_dof_count)
        r_vec_solved[free_dof_indices] = r_free_solved
        r_vec_solved[dependent_dof_indices] = x_d1.dot(r_free_solved) + x_d2.dot(r_fixed)
        r_vec_solved[fixed_dof_indices] = r_fixed

        print(f"Complete displacement vector: {r_vec_solved}")

        for element in self.truss.elements:
            element.set_local_deformations(r_vec_solved)
