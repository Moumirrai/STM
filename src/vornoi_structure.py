import math
from random import uniform
from time import perf_counter

import numpy as np
from scipy.spatial import Voronoi


def generateStructure(width, height, num_points, point_radius):
    # In rect with width and height, generate num_points random points.
    # Start with one random point, then for each next check if it is at least
    # 2 * point_radius away from all existing points; if not, generate a new point.
    # Limit attempts to 1000; if it fails to generate after that, return what we have.
    expandedDomainPoints = []
    innerPoints = []
    
    tile_offsets = (
        (0.0, 0.0),
        (width, 0.0),
        (-width, 0.0),
        (0.0, height),
        (0.0, -height),
        (width, height),
        (-width, height),
        (width, -height),
        (-width, -height),
    )
    
    attempts = 0
    min_distance_squared = (2 * point_radius) ** 2
    while len(innerPoints) < num_points and attempts < 1000:
        x = uniform(0, width)
        y = uniform(0, height)
        new_point = (x, y)
        if all(
            (x - p[0]) ** 2 + (y - p[1]) ** 2 >= min_distance_squared for p in expandedDomainPoints
        ):
            innerPoints.append(new_point)
            
            for dx, dy in tile_offsets:
                nx = x + dx
                ny = y + dy
                if nx >= 0 - 2*point_radius and nx <= width + 2*point_radius and ny >= 0 - 2*point_radius and ny <= height + 2*point_radius:
                    expandedDomainPoints.append((nx, ny))
            
            attempts = 0  # reset attempts after a successful addition
        else:
            attempts += 1

    

    tiled_points = []
    tiled_meta = []
    for base_index, (x, y) in enumerate(innerPoints):
        for dx, dy in tile_offsets:
            tiled_points.append((x + dx, y + dy))
            tiled_meta.append((base_index, dx, dy))

    vor = Voronoi(np.array(tiled_points))

    edges = []
    periodic_edges = []
    seen_edges = set()
    seen_periodic = set()
    print(vor.ridge_vertices[0])
    print(vor.ridge_points[0:5])
    print(vor.vertices[0:5])
    for (p1, p2), ridge_vertices in zip(vor.ridge_points, vor.ridge_vertices):
        if -1 in ridge_vertices:
            continue
        v0, v1 = vor.vertices[ridge_vertices]
        ridge_length = math.hypot(v1[0] - v0[0], v1[1] - v0[1])

        b1, dx1, dy1 = tiled_meta[p1]
        b2, dx2, dy2 = tiled_meta[p2]

        is_base_1 = dx1 == 0.0 and dy1 == 0.0
        is_base_2 = dx2 == 0.0 and dy2 == 0.0
        if not (is_base_1 or is_base_2):
            continue

        if is_base_1 and is_base_2:
            key = (min(b1, b2), max(b1, b2))
            if key not in seen_edges:
                edges.append((b1, b2, ridge_length))
                seen_edges.add(key)
            continue

        if is_base_1:
            base_idx = b1
            other_idx = b2
            shift = (dx2 - dx1, dy2 - dy1)
        else:
            base_idx = b2
            other_idx = b1
            shift = (dx1 - dx2, dy1 - dy2)

        key = (base_idx, other_idx, shift[0], shift[1])
        if key not in seen_periodic:
            periodic_edges.append((base_idx, other_idx, shift[0], shift[1], ridge_length))
            seen_periodic.add(key)

    return (innerPoints, edges, periodic_edges, expandedDomainPoints)

# Call function and plot with matplotlib.
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Rectangle

    width = 10
    height = 5
    num_points = 1000
    point_radius = 0.2
    start = perf_counter()
    points = generateStructure(width, height, num_points, point_radius)
    elapsed = perf_counter() - start
    innerPoints, edges, periodic_edges, expandedDomainPoints = points
    x, y = zip(*innerPoints)
    print(
        f"Generated {len(innerPoints)} inner points, {len(edges)} edges, "
        f"{len(periodic_edges)} periodic edges in {elapsed:.2f} seconds."
    )


    ax = plt.gca()
    ax.add_patch(Rectangle((0, 0), width, height, fill=False, linewidth=1))
    plt.scatter(x, y)
    for p1, p2, ridge_length in edges:
        x1, y1 = innerPoints[p1]
        x2, y2 = innerPoints[p2]
        ax.plot([x1, x2], [y1, y2], linewidth=max(0.5, ridge_length * 0.2))
    for p1, p2, dx, dy, ridge_length in periodic_edges:
        x1, y1 = innerPoints[p1]
        x2, y2 = innerPoints[p2]
        ax.plot(
            [x1, x2 + dx],
            [y1, y2 + dy],
            linewidth=max(0.5, ridge_length),
            linestyle="--",
        )
    for px, py in expandedDomainPoints:
        ax.add_patch(Circle((px, py), point_radius, fill=False, linewidth=1))
    plt.xlim(-width, 2 * width)
    plt.ylim(-height, 2 * height)
    ax.set_aspect("equal", adjustable="box")
    plt.title("Random Points with Minimum Distance")
    plt.show()