import math
from random import uniform
from time import perf_counter

import numpy as np
from scipy.spatial import Voronoi

from structure_parser import DependencyDefinition, EigenstrainDefinition, ElementDefinition, MasterDefinition, NodeDefinition, StructureDefinition


def _periodic_dist_sq(x1, y1, x2, y2, width, height, a = 1.0, b = 1.0):
    dx = min(abs(x1 - x2), width - abs(x1 - x2))
    dy = min(abs(y1 - y2), height - abs(y1 - y2))
    return ((dx * dx)/a) + ((dy * dy)/b)


def generateStructure(width, height, num_points, point_radius, a = 1.0, b = 1.0):
    # --- 1. Place points with minimum periodic separation (hard-disk packing) ---
    base_points = []
    min_dist_sq = (2 * point_radius) ** 2
    attempts = 0

    if min_dist_sq == 0.0:
        for _ in range(num_points):
            base_points.append((uniform(0, width), uniform(0, height)))
    else:
        min_dist = 2 * point_radius
        cell_size = min_dist
        nx = max(1, int(width / cell_size))
        ny = max(1, int(height / cell_size))
        grid = {}

        def _cell_index(x, y):
            # Clamp to avoid right/top edge mapping outside the last cell.
            ix = min(int(x / cell_size), nx - 1)
            iy = min(int(y / cell_size), ny - 1)
            return ix, iy

        while len(base_points) < num_points and attempts < 1000:
            x, y = uniform(0, width), uniform(0, height)
            ix, iy = _cell_index(x, y)

            valid = True
            for dix in (-1, 0, 1):
                for diy in (-1, 0, 1):
                    neighbor_key = ((ix + dix) % nx, (iy + diy) % ny)
                    for point_index in grid.get(neighbor_key, []):
                        px, py = base_points[point_index]
                        if _periodic_dist_sq(x, y, px, py, width, height, a, b) < min_dist_sq:
                            valid = False
                            break
                    if not valid:
                        break
                if not valid:
                    break

            if valid:
                point_index = len(base_points)
                base_points.append((x, y))
                grid.setdefault((ix, iy), []).append(point_index)
                attempts = 0
            else:
                attempts += 1

    # --- 2. Tile into 3x3 grid so every boundary cell is fully closed ---
    # Each base point is copied to the 8 surrounding tiles plus the base tile itself.
    # Every copy is tagged with its original index and the tile offset applied.
    tile_offsets = [
        (-width, -height), (0.0, -height), (width, -height),
        (-width,  0.0),    (0.0,  0.0),    (width,  0.0),
        (-width,  height), (0.0,  height), (width,  height),
    ]

    tiled_coords = []
    tiled_meta = []  # (original_index, tile_offset_x, tile_offset_y)
    for original_index, (x, y) in enumerate(base_points):
        for tile_offset_x, tile_offset_y in tile_offsets:
            tiled_coords.append((x + tile_offset_x, y + tile_offset_y))
            tiled_meta.append((original_index, tile_offset_x, tile_offset_y))

    vor = Voronoi(np.array(tiled_coords))

    # --- 3. Classify ridges and build output lists ---
    # A Voronoi ridge is the shared boundary between two seed-point cells.
    # We keep a ridge only when at least one seed is in the base tile (offset 0,0):
    #   both seeds in base tile → interior member between two base points
    #   one seed is a ghost     → periodic member; a ghost point is added to the list
    #
    # Ghost points represent where the neighbour lives when you look across the
    # boundary.  They are added to `points` and linked back to their original
    # base point via a dependency entry.

    points = [[x, y] for x, y in base_points]   # grows as ghost points are appended
    members = []        # [index_a, index_b, length]
    dependencies = []   # [master_index, ghost_index]  — ghost stitches to master
    seen = set()

    for (tiled_p, tiled_q), (ridge_vert_a, ridge_vert_b) in zip(vor.ridge_points, vor.ridge_vertices):
        if ridge_vert_a == -1 or ridge_vert_b == -1:   # open (infinite) ridge — skip
            continue

        origin_p, tile_dx_p, tile_dy_p = tiled_meta[tiled_p]
        origin_q, tile_dx_q, tile_dy_q = tiled_meta[tiled_q]

        is_base_p = tile_dx_p == 0.0 and tile_dy_p == 0.0
        is_base_q = tile_dx_q == 0.0 and tile_dy_q == 0.0

        if not (is_base_p or is_base_q):
            continue    # both seeds are ghosts — discard entirely

        ridge_length = math.dist(vor.vertices[ridge_vert_a], vor.vertices[ridge_vert_b])

        if is_base_p and is_base_q:
            # Interior member: both endpoints are base points
            key = (min(origin_p, origin_q), max(origin_p, origin_q))
            if key not in seen:
                members.append([origin_p, origin_q, ridge_length])
                seen.add(key)

        else:
            # Periodic member: one seed is a ghost copy from a neighbouring tile.
            # crossing_shift is the tile offset of the ghost seed — it tells us
            # which direction the boundary is crossed.
            if is_base_p:
                base_origin   = origin_p
                other_origin  = origin_q
                crossing_shift = (tile_dx_q, tile_dy_q)
            else:
                base_origin   = origin_q
                other_origin  = origin_p
                crossing_shift = (tile_dx_p, tile_dy_p)

            # The same physical crossing appears twice (once from each side of the
            # boundary) so check both forward and reverse keys before inserting.
            forward_key = (base_origin,  other_origin,  crossing_shift)
            reverse_key = (other_origin, base_origin,  (-crossing_shift[0], -crossing_shift[1]))

            if forward_key not in seen and reverse_key not in seen:
                ghost_x = base_points[other_origin][0] + crossing_shift[0]
                ghost_y = base_points[other_origin][1] + crossing_shift[1]
                ghost_index = len(points)
                points.append([ghost_x, ghost_y])
                members.append([base_origin, ghost_index, ridge_length])
                dependencies.append([other_origin, ghost_index])
                seen.add(forward_key)

    return points, members, dependencies

# Call function and plot with matplotlib.
if __name__ == "__main__REMOVE":
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Rectangle

    width = 10
    height = 5
    num_points = 25
    point_radius = 0.6
    start = perf_counter()
    points, members, dependencies = generateStructure(width, height, num_points, point_radius)
    elapsed = perf_counter() - start

    n_base = len(points) - len(dependencies)   # first n_base entries are base points

    print(f"Generated {n_base} base points + {len(dependencies)} ghost points, "
          f"{len(members)} members, {len(dependencies)} dependencies "
          f"in {elapsed:.2f} seconds.\n")

    print("Members [index_a, index_b, length]:")
    for m in members:
        print(f"  {m}")
    print("\nDependencies [master, ghost]:")
    for d in dependencies:
        print(f"  {d}")

    fig, ax = plt.subplots()
    ax.add_patch(Rectangle((0, 0), width, height, fill=False, linewidth=2, edgecolor="black"))

    # --- members: blue solid for interior, red dashed for periodic (ghost endpoint) ---
    for index_a, index_b, length in members:
        x1, y1 = points[index_a]
        x2, y2 = points[index_b]
        is_periodic = index_a >= n_base or index_b >= n_base
        style = dict(color="red", linestyle="--", linewidth=1) if is_periodic \
                else dict(color="blue", linestyle="-",  linewidth=1)
        ax.plot([x1, x2], [y1, y2], zorder=2)
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        label_color = "darkred" if is_periodic else "blue"
        ax.text(mx, my, f"{index_a}-{index_b}", fontsize=5, color=label_color,
                ha="center", va="center", zorder=7,
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.6))

    # --- base points: black dot + index label + exclusion-radius circle ---
    for index, (px, py) in enumerate(points[:n_base]):
        ax.plot(px, py, "ko", markersize=4, zorder=5)
        ax.text(px + 0.05, py + 0.05, str(index), fontsize=7, color="black", zorder=6)
        ax.add_patch(Circle((px, py), point_radius, fill=False, linewidth=0.5,
                             linestyle=":", edgecolor="gray"))

    # --- ghost points: red cross + index label + dependency annotation ---
    for master_index, ghost_index in dependencies:
        gx, gy = points[ghost_index]
        ax.plot(gx, gy, "rx", markersize=7, markeredgewidth=1.5, zorder=5)
        ax.text(gx + 0.05, gy + 0.05, str(ghost_index), fontsize=7, color="darkred", zorder=6)
        # Annotate which base point this ghost stitches to
        ax.text(gx + 0.05, gy - 0.25, f"→{master_index}", fontsize=5, color="darkred",
                zorder=6)

    ax.set_xlim(-width * 0.15, width * 1.15)
    ax.set_ylim(-height * 0.15, height * 1.15)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"{n_base} base pts · {len(dependencies)} ghosts · "
                 f"{len(members)} members · {len(dependencies)} deps")
    plt.tight_layout()
    plt.show()
    
def generate_voronoi_structure(width: float, height: float, num_points: int, point_radius: float, a: float = 1.0, b: float = 1.0) -> StructureDefinition:
    default_E = 30e9
    points, members, dependencies = generateStructure(width, height, num_points, point_radius, a, b)
    
    nodes = [NodeDefinition(dx=x, dy=y) for x, y in points]
    
    nodes[0].constraints = "xy"  # Fix the first node to prevent rigid body motion
    
    elements = []
    for index_a, index_b, length in members:
        # Cross-section proportional to Voronoi ridge length; E is set globally via defaultYoungsModulus
        elements.append(ElementDefinition(
            starting_node=index_a,
            ending_node=index_b,
            A=length,
        ))
        
    dependencies = [DependencyDefinition(
        node=ghost_index,
        masters=[MasterDefinition(node=master_index, direction="x", factor=1.0), MasterDefinition(node=master_index, direction="y", factor=1.0)]
    ) for master_index, ghost_index in dependencies]
    
    
        
    return StructureDefinition(
        nodes=nodes,
        elements=elements,
        dependencies=dependencies,
        eigenstrain=EigenstrainDefinition(x=1.0, y=0.0, angle=1.0),
        defaultYoungsModulus=default_E,
        volume=width * height,
    )

def aniotropise(structure: StructureDefinition, xE: float, yE: float) -> StructureDefinition:
    for element in structure.elements:
        node_a = structure.nodes[element.starting_node]
        node_b = structure.nodes[element.ending_node]
        dx = node_b.dx - node_a.dx
        dy = node_b.dy - node_a.dy
        angle = math.atan2(dy, dx)
        cos_angle = math.cos(angle)
        sin_angle = math.sin(angle)
        effective_E = 1 / ((cos_angle ** 2) / xE + (sin_angle ** 2) / yE)
        element.E = effective_E
    return structure