import numpy as np

import matplotlib
matplotlib.use('qtagg')

from models import TrussData
from termcolor import colored
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


def _select_tile_vectors(truss: TrussData, deformed_positions: dict[int, tuple[float, float]]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Determine two basis vectors for tiling using dependency pairs when possible.
    Returns (before_v1, before_v2, after_v1, after_v2).
    """
    pair_vectors: list[tuple[int, int, np.ndarray, np.ndarray]] = []

    for node in truss.nodes:
        if node.dependency is None or not node.dependency.masters:
            continue

        # A dependent node may have multiple DOF masters but often maps to the same master node.
        unique_master_indices: list[int] = []
        for master in node.dependency.masters:
            if master.nodeIndex not in unique_master_indices:
                unique_master_indices.append(master.nodeIndex)

        for master_index in unique_master_indices:
            master_node = truss.nodes[master_index]
            before_vec = np.array([node.dx - master_node.dx, node.dy - master_node.dy], dtype=float)
            if np.linalg.norm(before_vec) == 0.0:
                continue

            dep_def = np.array(deformed_positions[node.index], dtype=float)
            master_def = np.array(deformed_positions[master_index], dtype=float)
            after_vec = dep_def - master_def
            pair_vectors.append((master_index, node.index, before_vec, after_vec))

    if pair_vectors:
        print(colored("Dependency vectors (master -> dependent):", "magenta"))
        for master_idx, dep_idx, before_vec, after_vec in pair_vectors:
            print(
                colored(
                    f"  {master_idx} -> {dep_idx}: before {tuple(before_vec)}, after {tuple(after_vec)}",
                    "magenta",
                )
            )

    # Prefer dependency-derived vectors and pick a non-collinear pair.
    if len(pair_vectors) >= 2:
        pair_vectors.sort(key=lambda item: np.linalg.norm(item[2]), reverse=False)
        _, _, before_v1, after_v1 = pair_vectors[0]

        chosen_second: tuple[np.ndarray, np.ndarray] | None = None
        for _, _, candidate_before, candidate_after in pair_vectors[1:]:
            denom = np.linalg.norm(before_v1) * np.linalg.norm(candidate_before)
            if denom == 0.0:
                continue
            cross = abs(before_v1[0] * candidate_before[1] - before_v1[1] * candidate_before[0])
            if cross / denom > 1e-3:
                chosen_second = (candidate_before, candidate_after)
                break

        if chosen_second is not None:
            before_v2, after_v2 = chosen_second
            return before_v1, before_v2, after_v1, after_v2

    # Fallback to axis-aligned cell extents if dependency vectors are missing/collinear.
    x_coords = [node.dx for node in truss.nodes]
    y_coords = [node.dy for node in truss.nodes]
    dx = max(x_coords) - min(x_coords) if x_coords else 1.0
    dy = max(y_coords) - min(y_coords) if y_coords else 1.0
    if dx == 0.0:
        dx = 1.0
    if dy == 0.0:
        dy = 1.0

    def_x = [coord[0] for coord in deformed_positions.values()]
    def_y = [coord[1] for coord in deformed_positions.values()]
    dfx = max(def_x) - min(def_x) if def_x else dx
    dfy = max(def_y) - min(def_y) if def_y else dy
    if dfx == 0.0:
        dfx = dx
    if dfy == 0.0:
        dfy = dy

    before_v1 = np.array([dx, 0.0], dtype=float)
    before_v2 = np.array([0.0, dy], dtype=float)
    after_v1 = np.array([dfx, 0.0], dtype=float)
    after_v2 = np.array([0.0, dfy], dtype=float)
    return before_v1, before_v2, after_v1, after_v2

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


def plot_deformed_structure(truss: TrussData, scale: float = 1, tiled: bool = False, tile_opacity: float = 0.25, original: bool = False):
    # Compute deformed positions
    deformed_positions: dict[int, tuple[float, float]] = {}
    for node in truss.nodes:
        deltax = node.local_deformations[0] * scale if node.local_deformations is not None else 0.0
        deltay = node.local_deformations[1] * scale if node.local_deformations is not None else 0.0
        x = node.dx + deltax
        y = node.dy + deltay
        deformed_positions[node.index] = (x, y)
        #print(colored(f"Node {node.index}: Original ({node.dx}, {node.dy}), Deformed by ({deltax}, {deltay})", "cyan"))

    # Get axial forces for coloring
    forces = [element.axial_force() for element in truss.elements]
    
    #print(colored(f"Axial forces for coloring: {forces}", "yellow"))

    if not forces:
        print(colored("No elements found. Nothing to plot.", "red"))
        return

    # Normalize forces for colormap
    min_force = min(forces)
    max_force = max(forces)
    if min_force == max_force:
        # Keep a tiny range so colormap normalization is well-defined.
        min_force -= 1e-12
        max_force += 1e-12
    norm = Normalize(vmin=min_force, vmax=max_force)
    cmap = plt.get_cmap('viridis')

    fig, ax = plt.subplots()

    if tiled:
        before_v1, before_v2, after_v1, after_v2 = _select_tile_vectors(truss, deformed_positions)
        print(
            colored(
                f"Tile vectors - before: {tuple(before_v1)}, {tuple(before_v2)} | after: {tuple(after_v1)}, {tuple(after_v2)}",
                "green",
            )
        )

        tile_offsets: list[tuple[float, float, float]] = []
        for iy in (-1, 0, 1):
            for ix in (-1, 0, 1):
                shift = ix * after_v1 + iy * after_v2
                alpha = 1.0 if ix == 0 and iy == 0 else tile_opacity
                tile_offsets.append((shift[0], shift[1], alpha))
    else:
        tile_offsets = [(0.0, 0.0, 1.0)]

    # Plot undeformed structure in background
    if original:
        for element in truss.elements:
            n1, n2 = element.nodes
            x1, y1 = n1.dx, n1.dy
            x2, y2 = n2.dx, n2.dy
            ax.plot([x1, x2], [y1, y2], color='black', alpha=0.35, linewidth=2)

    # Plot each element with color based on axial force
    for xoff, yoff, alpha in tile_offsets:
        for element, force in zip(truss.elements, forces):
            n1, n2 = element.nodes
            x1, y1 = deformed_positions[n1.index]
            x2, y2 = deformed_positions[n2.index]
            color = cmap(norm(force))
            ax.plot([x1 + xoff, x2 + xoff], [y1 + yoff, y2 + yoff], color=color, alpha=alpha, linewidth=2)

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Axial Force')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Deformed Structure' + (' (3x3 Tiled)' if tiled else ''))
    ax.axis('equal')
    plt.show()
    
def plot_deformed_structure_black(truss: TrussData, scale: float = 1, tiled: bool = False, tile_opacity: float = 0.25, original: bool = False):
    # Compute deformed positions
    deformed_positions: dict[int, tuple[float, float]] = {}
    for node in truss.nodes:
        deltax = node.local_deformations[0] * scale if node.local_deformations is not None else 0.0
        deltay = node.local_deformations[1] * scale if node.local_deformations is not None else 0.0
        x = node.dx + deltax
        y = node.dy + deltay
        deformed_positions[node.index] = (x, y)
        #print(colored(f"Node {node.index}: Original ({node.dx}, {node.dy}), Deformed by ({deltax}, {deltay})", "cyan"))

    # Get axial forces for coloring
    forces = [element.axial_force() for element in truss.elements]
    
    #print(colored(f"Axial forces for coloring: {forces}", "yellow"))

    if not forces:
        print(colored("No elements found. Nothing to plot.", "red"))
        return

    # Normalize forces for colormap
    min_force = min(forces)
    max_force = max(forces)
    if min_force == max_force:
        # Keep a tiny range so colormap normalization is well-defined.
        min_force -= 1e-12
        max_force += 1e-12
    norm = Normalize(vmin=min_force, vmax=max_force)
    cmap = plt.get_cmap('viridis')

    fig, ax = plt.subplots()

    if tiled:
        before_v1, before_v2, after_v1, after_v2 = _select_tile_vectors(truss, deformed_positions)
        print(
            colored(
                f"Tile vectors - before: {tuple(before_v1)}, {tuple(before_v2)} | after: {tuple(after_v1)}, {tuple(after_v2)}",
                "green",
            )
        )

        tile_offsets: list[tuple[float, float, float]] = []
        for iy in (-1, 0, 1):
            for ix in (-1, 0, 1):
                shift = ix * after_v1 + iy * after_v2
                alpha = 1.0 if ix == 0 and iy == 0 else tile_opacity
                tile_offsets.append((shift[0], shift[1], alpha))
    else:
        tile_offsets = [(0.0, 0.0, 1.0)]

    # Plot undeformed structure in background
    if original:
        for element in truss.elements:
            n1, n2 = element.nodes
            x1, y1 = n1.dx, n1.dy
            x2, y2 = n2.dx, n2.dy
            ax.plot([x1, x2], [y1, y2], color='black', alpha=0.35, linewidth=2)

    # Plot each element with color based on axial force
    for xoff, yoff, alpha in tile_offsets:
        for element, force in zip(truss.elements, forces):
            n1, n2 = element.nodes
            x1, y1 = deformed_positions[n1.index]
            x2, y2 = deformed_positions[n2.index]
            color = cmap(norm(force))
            ax.plot([x1 + xoff, x2 + xoff], [y1 + yoff, y2 + yoff], color='black', alpha=alpha, linewidth=2)

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Axial Force')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Deformed Structure' + (' (3x3 Tiled)' if tiled else ''))
    ax.axis('equal')
    plt.show()