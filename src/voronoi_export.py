import math
import argparse
from pathlib import Path
from random import Random, uniform
from time import perf_counter
from typing import Any, Callable

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import numpy as np
from scipy.spatial import Delaunay, Voronoi

from vornoi_structure import _periodic_dist_sq


def _generate_base_points(
    width: float,
    height: float,
    num_points: int,
    point_radius: float,
    a: float = 1.0,
    b: float = 1.0,
    uniform_fn: Callable[[float, float], float] = uniform,
) -> list[tuple[float, float]]:
    base_points: list[tuple[float, float]] = []
    min_dist_sq = (2 * point_radius) ** 2
    attempts = 0

    if min_dist_sq == 0.0:
        for _ in range(num_points):
            base_points.append((float(uniform_fn(0, width)), float(uniform_fn(0, height))))
        return base_points

    min_dist = 2 * point_radius
    cell_size = min_dist
    nx = max(1, int(width / cell_size))
    ny = max(1, int(height / cell_size))
    grid: dict[tuple[int, int], list[int]] = {}

    def _cell_index(x: float, y: float) -> tuple[int, int]:
        ix = min(int(x / cell_size), nx - 1)
        iy = min(int(y / cell_size), ny - 1)
        return ix, iy

    while len(base_points) < num_points and attempts < 1000:
        x = float(uniform_fn(0, width))
        y = float(uniform_fn(0, height))
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

    return base_points


def _tile_points(
    base_points: list[tuple[float, float]],
    width: float,
    height: float,
) -> tuple[list[tuple[float, float]], list[tuple[int, float, float]]]:
    tile_offsets = [
        (-width, -height), (0.0, -height), (width, -height),
        (-width, 0.0), (0.0, 0.0), (width, 0.0),
        (-width, height), (0.0, height), (width, height),
    ]

    tiled_coords: list[tuple[float, float]] = []
    tiled_meta: list[tuple[int, float, float]] = []
    for original_index, (x, y) in enumerate(base_points):
        for tile_offset_x, tile_offset_y in tile_offsets:
            tiled_coords.append((x + tile_offset_x, y + tile_offset_y))
            tiled_meta.append((original_index, tile_offset_x, tile_offset_y))

    return tiled_coords, tiled_meta


def _build_voronoi_geometry(
    width: float,
    height: float,
    num_points: int,
    point_radius: float,
    a: float = 1.0,
    b: float = 1.0,
    uniform_fn: Callable[[float, float], float] = uniform,
) -> tuple[list[tuple[float, float]], list[tuple[int, int, float]], list[tuple[int, int, float, float, float]], Voronoi]:
    base_points = _generate_base_points(width, height, num_points, point_radius, a, b, uniform_fn)
    tiled_coords, tiled_meta = _tile_points(base_points, width, height)
    vor = Voronoi(np.array(tiled_coords))

    members: list[tuple[int, int, float]] = []
    periodic_members: list[tuple[int, int, float, float, float]] = []
    seen: set[tuple[Any, ...]] = set()

    for (tiled_p, tiled_q), (ridge_vert_a, ridge_vert_b) in zip(vor.ridge_points, vor.ridge_vertices):
        if ridge_vert_a == -1 or ridge_vert_b == -1:
            continue

        origin_p, tile_dx_p, tile_dy_p = tiled_meta[tiled_p]
        origin_q, tile_dx_q, tile_dy_q = tiled_meta[tiled_q]

        is_base_p = tile_dx_p == 0.0 and tile_dy_p == 0.0
        is_base_q = tile_dx_q == 0.0 and tile_dy_q == 0.0

        if not (is_base_p or is_base_q):
            continue

        ridge_length = float(math.dist(vor.vertices[ridge_vert_a], vor.vertices[ridge_vert_b]))

        if is_base_p and is_base_q:
            key = (min(origin_p, origin_q), max(origin_p, origin_q))
            if key not in seen:
                members.append((origin_p, origin_q, ridge_length))
                seen.add(key)
        else:
            if is_base_p:
                base_origin = origin_p
                other_origin = origin_q
                crossing_shift = (tile_dx_q, tile_dy_q)
            else:
                base_origin = origin_q
                other_origin = origin_p
                crossing_shift = (tile_dx_p, tile_dy_p)

            forward_key = (base_origin, other_origin, crossing_shift)
            reverse_key = (other_origin, base_origin, (-crossing_shift[0], -crossing_shift[1]))

            if forward_key not in seen and reverse_key not in seen:
                periodic_members.append(
                    (
                        base_origin,
                        other_origin,
                        float(crossing_shift[0]),
                        float(crossing_shift[1]),
                        ridge_length,
                    )
                )
                seen.add(forward_key)

    return base_points, members, periodic_members, vor


def _plot_point_view(
    output_path: Path,
    base_points: list[tuple[float, float]],
    width: float,
    height: float,
    point_radius: float,
) -> Path:
    figure, axis = plt.subplots(figsize=(10, 10), dpi=100)
    axis.add_patch(Rectangle((0, 0), width, height, fill=False, linewidth=2, edgecolor="black"))

    for x, y in base_points:
        axis.plot(x, y, "ko", markersize=8, zorder=5)
        axis.add_patch(
            Circle((x, y), point_radius, fill=False, linewidth=2, linestyle="--", edgecolor="gray")
        )

    axis.set_xlim(-width * 0.15, width * 1.15)
    axis.set_ylim(-height * 0.15, height * 1.15)
    axis.set_aspect("equal", adjustable="box")
    axis.set_axis_off()
    figure.subplots_adjust(left=0, right=1, bottom=0, top=1)
    figure.savefig(output_path, dpi=100, pad_inches=0)
    plt.close(figure)
    return output_path


def _plot_geometry_view(
    output_path: Path,
    base_points: list[tuple[float, float]],
    members: list[tuple[int, int, float]],
    periodic_members: list[tuple[int, int, float, float, float]],
    vor: Voronoi,
    width: float,
    height: float,
) -> Path:
    figure, axis = plt.subplots(figsize=(10, 10), dpi=100)
    axis.add_patch(Rectangle((0, 0), width, height, fill=False, linewidth=2, edgecolor="black"))

    delaunay = Delaunay(np.array(base_points))
    for simplex in delaunay.simplices:
        triangle_indices = list(simplex) + [int(simplex[0])]
        triangle = np.array([base_points[index] for index in triangle_indices])
        #axis.plot(triangle[:, 0], triangle[:, 1], color="0.75", linewidth=0.8, zorder=1)

    for vertex_a, vertex_b in vor.ridge_vertices:
        if vertex_a == -1 or vertex_b == -1:
            continue
        segment = vor.vertices[[vertex_a, vertex_b]]
        axis.plot(segment[:, 0], segment[:, 1], color="tab:red", linewidth=1.0, zorder=2)

    for index_a, index_b, _length in members:
        x1, y1 = base_points[index_a]
        x2, y2 = base_points[index_b]
        axis.plot([x1, x2], [y1, y2], color="black", linewidth=3, zorder=3)

    for index_a, index_b, shift_x, shift_y, _length in periodic_members:
        x1, y1 = base_points[index_a]
        x2, y2 = base_points[index_b]
        axis.plot(
            [x1, x2 + shift_x],
            [y1, y2 + shift_y],
            color="black",
            linestyle="--",
            linewidth=3,
            zorder=3,
        )

    for x, y in base_points:
        axis.plot(x, y, "ko", markersize=8, zorder=4)

    axis.set_xlim(-width * 0.15, width * 1.15)
    axis.set_ylim(-height * 0.15, height * 1.15)
    axis.set_aspect("equal", adjustable="box")
    axis.set_axis_off()
    figure.subplots_adjust(left=0, right=1, bottom=0, top=1)
    figure.savefig(output_path, dpi=100, pad_inches=0)
    plt.close(figure)
    return output_path


def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Render periodic Voronoi geometry as matplotlib figures.")
    parser.add_argument("width", type=float, help="Domain width.")
    parser.add_argument("height", type=float, help="Domain height.")
    parser.add_argument("num_points", type=int, help="Number of random seed points to place.")
    parser.add_argument("point_radius", type=float, help="Hard-disk exclusion radius for points.")
    parser.add_argument("--a", type=float, default=1.0, help="Anisotropy factor in x.")
    parser.add_argument("--b", type=float, default=1.0, help="Anisotropy factor in y.")
    parser.add_argument("--seed", type=int, default=None, help="Optional random seed for reproducible output.")
    parser.add_argument(
        "--points-output",
        type=Path,
        default=Path("voronoi_points.png"),
        help="Output image path for the seed-point view.",
    )
    parser.add_argument(
        "--geometry-output",
        type=Path,
        default=Path("voronoi_geometry.png"),
        help="Output image path for the triangulation/ridge view.",
    )
    return parser


def main() -> None:
    parser = _build_argument_parser()
    args = parser.parse_args()

    if args.width <= 0 or args.height <= 0:
        raise ValueError("width and height must be positive")
    if args.num_points <= 0:
        raise ValueError("num_points must be positive")
    if args.point_radius < 0:
        raise ValueError("point_radius must be non-negative")

    start = perf_counter()
    if args.seed is None:
        uniform_fn = uniform
    else:
        rng = Random(args.seed)
        uniform_fn = rng.uniform

    base_points, members, periodic_members, vor = _build_voronoi_geometry(
        width=args.width,
        height=args.height,
        num_points=args.num_points,
        point_radius=args.point_radius,
        a=args.a,
        b=args.b,
        uniform_fn=uniform_fn,
    )

    points_output = _plot_point_view(
        output_path=args.points_output,
        base_points=base_points,
        width=args.width,
        height=args.height,
        point_radius=args.point_radius,
    )
    geometry_output = _plot_geometry_view(
        output_path=args.geometry_output,
        base_points=base_points,
        members=members,
        periodic_members=periodic_members,
        vor=vor,
        width=args.width,
        height=args.height,
    )
    elapsed = perf_counter() - start
    print(f"Saved {points_output} and {geometry_output} in {elapsed:.2f} s")


if __name__ == "__main__":
    main()