import numpy as np

from models import TrussData
import pyvista as pv

def export_vtk(truss: TrussData):
    # convert nodes and deformations to Vec3
    points = np.array([[node.dx, node.dy, 0.0] for node in truss.nodes])
    displacements = np.array([
        [node.local_deformations[0] if node.local_deformations is not None else 0.0,
        node.local_deformations[1] if node.local_deformations is not None else 0.0,
        0.0] for node in truss.nodes
    ])

    forces = np.array([element.axial_force() for element in truss.elements])

    # create lines from elements, 2 specifies number of points per line
    lines = np.array([[2, element.nodes[0].index, element.nodes[1].index] for element in truss.elements]).flatten()

    poly = pv.PolyData()
    poly.points = points
    poly.lines = lines
    poly.point_data['displacement'] = displacements
    poly.cell_data['axial_force'] = forces

    poly.save("fixed.vtp", binary=True)