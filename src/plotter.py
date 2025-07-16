import matplotlib.pyplot as plt
import matplotlib.colors as matcolor
import numpy as np

from models import TrussData
import pyvista as pv


def render_truss_structure(truss: TrussData):
    for i, node in enumerate(truss.nodes):
        plt.plot(node.dx, node.dy, 'o', color='black')
        plt.text(node.dx, node.dy - 0.02, f"N{i}", fontsize=8, ha='right')

    for element in truss.elements:
        plt.plot([element.nodes[0].dx, element.nodes[1].dx], [element.nodes[0].dy, element.nodes[1].dy], 'k-')

    plt.xlabel('dx')
    plt.ylabel('dy')
    plt.title('Nodes and Elements')
    plt.grid(True)

    plt.gca().invert_yaxis()
    plt.show()


def visualise_axial_forces(truss: TrussData):
    axial_forces: list[float] = [0.0] * len(truss.elements)
    for i, element in enumerate(truss.elements):
        axial_forces[i] = element.axial_force()

    min_axial_force = min(axial_forces)
    max_axial_force = max(axial_forces)

    if min_axial_force == max_axial_force:
        norm = matcolor.Normalize(vmin=min_axial_force - 0.1,
                                  vmax=max_axial_force + 0.1)  # add a small range to avoid zero division
    else:
        norm = matcolor.Normalize(vmin=min_axial_force, vmax=max_axial_force)

    colormap = plt.cm.plasma
    sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
    sm.set_array([])

    fig, ax = plt.subplots()

    for i, element in enumerate(truss.elements):
        force = axial_forces[i]
        ax.plot(
            [element.nodes[0].dx, element.nodes[1].dx],
            [element.nodes[0].dy, element.nodes[1].dy],
            color=colormap(norm(force)),
            linewidth=2,
        )
        ax.text(
            (element.nodes[0].dx + element.nodes[1].dx) / 2,
            (element.nodes[0].dy + element.nodes[1].dy) / 2,
            f"{force:.2f} N",
            fontsize=8,
            ha='center',
        )

    for node in truss.nodes:
        ax.plot(node.dx, node.dy, 'o', color='black')

    colorbar = fig.colorbar(sm, ax=ax, shrink=0.7)
    colorbar.set_label('Axial Force (N)')
    ax.invert_yaxis()
    ax.set_title('Axial Forces in Elements')
    ax.grid(False)

    fig.tight_layout()

    plt.show()


def export_vtk(truss: TrussData):
    # convert nodes and deformations to Vec3
    points = np.array([[node.dx, node.dy, 0.0] for node in truss.nodes])
    displacements = np.array([[node.local_deformations[0], node.local_deformations[1], 0.0] for node in truss.nodes])

    forces = np.array([element.axial_force() for element in truss.elements])

    # create lines from elements, 2 specifies number of points per line
    lines = np.array([[2, element.nodes[0].index, element.nodes[1].index] for element in truss.elements]).flatten()

    poly = pv.PolyData()
    poly.points = points
    poly.lines = lines
    poly.point_data['displacement'] = displacements
    poly.cell_data['axial_force'] = forces

    poly.save("new.vtp", binary=True)
