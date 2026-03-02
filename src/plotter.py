import numpy as np

from models import TrussData
from termcolor import colored
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

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


def export_vtk_to_file(truss: TrussData, filename: str):
    # convert nodes and deformations to Vec3
    points = np.array([[node.dx, node.dy, 0.0] for node in truss.nodes])
    displacements = np.array([
        [node.local_deformations[0] if node.local_deformations is not None else 0.0,
        node.local_deformations[1] if node.local_deformations is not None else 0.0,
        0.0] for node in truss.nodes
    ])

    with open(filename, 'w') as f:
        f.write("#let points = (\n")
        for node in truss.nodes:
            f.write(f"    ({node.dx}, {node.dy}),\n")
        f.write(")\n")

        f.write("#let connections = (\n")
        for element in truss.elements:
            n1 = element.nodes[0].index
            n2 = element.nodes[1].index
            f.write(f"    (\"{n1}\", \"{n2}\"),\n")
        f.write(")\n")

        f.write('#let displacements = (\n')
        for dx, dy, _ in displacements:
            f.write(f"    ({dx}, {dy}),\n")
        f.write(")\n")

    forces = np.array([element.axial_force() for element in truss.elements])

    # create lines from elements, 2 specifies number of points per line
    lines = np.array([[2, element.nodes[0].index, element.nodes[1].index] for element in truss.elements]).flatten()

    # construct PolyData with points and lines to avoid assignment-type mismatch


def plot_deformed_structure(truss: TrussData):
    scale = 0.5
    # Compute deformed positions
    deformed_positions = []
    for node in truss.nodes:
        deltax = node.local_deformations[0] * scale if node.local_deformations is not None else 0.0
        deltay = node.local_deformations[1] * scale if node.local_deformations is not None else 0.0
        x = node.dx + deltax
        y = node.dy + deltay
        deformed_positions.append((x, y))
        print(colored(f"Node {node.index}: Original ({node.dx}, {node.dy}), Deformed by ({deltax}, {deltay})", "cyan"))

    # Get axial forces for coloring
    forces = [element.axial_force() for element in truss.elements]
    
    print(colored(f"Axial forces for coloring: {forces}", "yellow"))

    # Normalize forces for colormap
    norm = Normalize(vmin=min(forces), vmax=max(forces))
    cmap = plt.get_cmap('viridis')

    fig, ax = plt.subplots()

    # Plot undeformed structure in background
    for element in truss.elements:
        n1, n2 = element.nodes
        x1, y1 = n1.dx, n1.dy
        x2, y2 = n2.dx, n2.dy
        ax.plot([x1, x2], [y1, y2], color='black', alpha=0.5, linewidth=2)

    # Plot each element with color based on axial force
    for element, force in zip(truss.elements, forces):
        n1, n2 = element.nodes
        x1, y1 = deformed_positions[n1.index]
        x2, y2 = deformed_positions[n2.index]
        color = cmap(norm(force))
        ax.plot([x1, x2], [y1, y2], color=color, linewidth=2)

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Axial Force')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Deformed Structure')
    ax.axis('equal')
    plt.show()