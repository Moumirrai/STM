import numpy as np

from models import TrussData
import pyvista as pv
from termcolor import colored

def export_vtk(truss: TrussData):
    # convert nodes and deformations to Vec3
    points = np.array([[node.dx, node.dy, 0.0] for node in truss.nodes])
    displacements = np.array([
        [node.local_deformations[0] if node.local_deformations is not None else 0.0,
        node.local_deformations[1] if node.local_deformations is not None else 0.0,
        0.0] for node in truss.nodes
    ])


    print(colored("#let points = (","black", "on_light_blue"))
    for node in truss.nodes:
        print(colored(f"    ({node.dx}, {node.dy}),", "light_blue"))
    print(colored(")", "light_blue"))

    print(colored("#let connections = (","black", "on_light_yellow"))
    for element in truss.elements:
        n1 = element.nodes[0].index
        n2 = element.nodes[1].index
        print(colored(f"    (\"{n1}\", \"{n2}\"),", "light_yellow"))
    print(colored(")", "light_yellow"))

    print(colored('#let displacements = (', 'black', 'on_light_green'))
    for dx, dy, _ in displacements:
        print(colored(f"    ({dx}, {dy}),", "light_green"))
    print(colored(")", "light_green"))


    forces = np.array([element.axial_force() for element in truss.elements])

    # create lines from elements, 2 specifies number of points per line
    lines = np.array([[2, element.nodes[0].index, element.nodes[1].index] for element in truss.elements]).flatten()

    # construct PolyData with points and lines to avoid assignment-type mismatch
    poly = pv.PolyData(points, lines)
    poly.point_data['displacement'] = displacements
    poly.cell_data['axial_force'] = forces

    poly.save("output.vtp", binary=True)