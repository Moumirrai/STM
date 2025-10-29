import numpy as np
from models import TrussData
from scipy.sparse import lil_matrix, linalg, identity, bmat, csr_matrix

from scipy.sparse.linalg import spsolve

E = 210e6  # Pa
A = 0.01  # m^2


class TrussSolver:

    def __init__(
            self, truss: TrussData, young_modulus: float = 210e6, area: float = 0.01
    ):
        self.truss = truss # degrees of freedom IDs

    def solve(self):

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
        DOFIDS = free_dof_indices + dependent_dof_indices + fixed_dof_indices

        """
        X matrix structure:
          +-------+-------+
          │  X11  │  None │ <- Free DOFs
          +-------+-------+
          │  XD1  │  XD2  │ <- Dependent DOFs
          +-------+-------+
          │  None │  X22  │ <- Fixed DOFs
          +-------+-------+
             Free   Fixed
        """

        X11 = identity(len(free_dof_indices))
        X22 = identity(len(fixed_dof_indices))

        free_and_fixed_indices = free_dof_indices + fixed_dof_indices

        # create the xd1+xd2 matrix
        XD = lil_matrix((len(dependent_dof_indices), len(free_and_fixed_indices)))
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
                    XD[local_i, master_local_index] = master.factor

        XD = XD.tocsr()
        # now we can split the XD matrix into XD1 (free DOFs) and XD2 (fixed DOFs)
        XD1 = XD[:, :len(free_dof_indices)]
        XD2 = XD[:, len(free_dof_indices):]

        # and assemble the full X matrix
        x_mat = bmat([
            [X11, None],
            [XD1, XD2],
            [None, X22]
        ])
        x_mat = x_mat.tocsr()
        
        # initialize reduced displacement and force vectors with known lengths
        u_reduced = np.zeros(len(free_and_fixed_indices))
        f_vec = np.zeros(len(free_dof_indices) + len(dependent_dof_indices))
        a_dependant_vec = np.zeros(len(dependent_dof_indices))

        free_index_map = {dof: idx for idx, dof in enumerate(free_dof_indices)}
        fixed_index_map = {dof: idx for idx, dof in enumerate(fixed_dof_indices)}
        dependent_index_map = {dof: idx for idx, dof in enumerate(dependent_dof_indices)}

        for node_idx, node in enumerate(self.truss.nodes):
            base_dof = node_idx * 2
            loads = [node.load_x, node.load_y]
            deformations = [node.deformation_x, node.deformation_y]
            eigenstrains = node.eigenstrain

            for direction, (load, deformation, eigenstrain) in enumerate(zip(loads, deformations, eigenstrains)):
                global_dof = base_dof + direction

                if global_dof in free_index_map:
                    u_reduced[free_index_map[global_dof]] = deformation
                    f_vec[free_index_map[global_dof]] = load
                elif global_dof in dependent_index_map:
                    f_vec[len(free_dof_indices) + dependent_index_map[global_dof]] = load
                    a_dependant_vec[dependent_index_map[global_dof]] = eigenstrain
                elif global_dof in fixed_index_map:
                    u_reduced[len(free_dof_indices) + fixed_index_map[global_dof]] = deformation

        u_vec = x_mat.dot(u_reduced)

        # split u_vec into free and fixed parts
        u_fixed = u_vec[len(free_dof_indices) + len(dependent_dof_indices):]

        f_1 = f_vec[:len(free_dof_indices)] #TODO
        f_D = f_vec[len(free_dof_indices):len(free_dof_indices) + len(dependent_dof_indices)] #TODO

        # now we need to assamble the global stiffness matrix K
        # we will devide it into 3 parts:
        # K11 - free DOFs, K1D - dependent DOFs, K12 - fixed DOFs, only the horizontal part is needed

        raw_K_matrix = lil_matrix((total_dof_count, total_dof_count))

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

            # global DOF indices for this element
            dofs = element.getDOFs()

            # add to global stiffness matrix
            for i in range(4):
                for j in range(4):
                    raw_K_matrix[dofs[i], dofs[j]] += stiffness_matrix[i, j]

        raw_K_matrix = raw_K_matrix.tocsr()

        # Partition K matrix according to DOF ordering
        # convert index lists to numpy int arrays for indexing
        free_idx = np.array(free_dof_indices, dtype=int)
        dep_idx = np.array(dependent_dof_indices, dtype=int)
        fix_idx = np.array(fixed_dof_indices, dtype=int)

        K11 = raw_K_matrix[free_idx][:, free_idx]
        K12 = raw_K_matrix[free_idx][:, fix_idx]
        K1D = raw_K_matrix[free_idx][:, dep_idx]
        KD1 = raw_K_matrix[dep_idx][:, free_idx]
        KDD = raw_K_matrix[dep_idx][:, dep_idx]
        K21 = raw_K_matrix[fix_idx][:, free_idx]
        KD2 = raw_K_matrix[dep_idx][:, fix_idx]
        K2D = raw_K_matrix[fix_idx][:, dep_idx]

        ASSAMBLED_K = K11 + XD1.T @ KD1 + K1D @ XD1 + XD1.T @ KDD @ XD1

        ASSAMBLED_F = -1 * (f_1 + XD1.T @ f_D + (0.5 * u_fixed.T @ (XD2.T @ KD1 + K21 + KD2.T @ KDD @ KD1+K2D @ XD1)).T + 0.5 * (K1D @ XD2 + XD1.T @ KDD @ XD2 + K12 + XD1.T @ KD2) @ u_fixed + (K1D + XD1.T @ KDD) @ a_dependant_vec)

        u_free_solved = spsolve(ASSAMBLED_K, ASSAMBLED_F)

        print(f"Solved free displacements: {u_free_solved}")

        # Update the full displacement vector
        u_vec_solved = np.zeros(total_dof_count)
        u_vec_solved[free_dof_indices] = u_free_solved
        u_vec_solved[dependent_dof_indices] = XD1.dot(u_free_solved) + XD2.dot(u_fixed) + a_dependant_vec
        u_vec_solved[fixed_dof_indices] = np.array(u_fixed).flatten()

        print(f"Complete displacement vector: {u_vec_solved}")

        for element in self.truss.elements:
            element.set_local_deformations(u_vec_solved)
