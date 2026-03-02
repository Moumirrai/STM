# Periodic Voronoi Structure Generator

## Overview

`vornoi_structure.py` generates a random 2-D point set inside a rectangular
domain and from it builds a **Voronoi tessellation** that tiles perfectly under
periodic (wrap-around) boundary conditions.  The output is a graph whose nodes
are the random points and whose edges are the Voronoi ridges, split into
*interior edges* (both endpoints inside the base cell) and *periodic edges*
(the ridge crosses the cell boundary).

---

## Concepts

### The 2-D Torus

A *periodic* or *toroidal* domain of width $W$ and height $H$ identifies the
left edge with the right edge and the top edge with the bottom edge.  Any
structure built on such a domain tiles the plane seamlessly: placing copies at
offsets $(nW, mH)$ for integer $n, m$ joins up without gaps or overlaps.

### Voronoi Tessellation

Given a set of seed points $\{s_i\}$, the Voronoi cell of $s_i$ is the region
of the plane closer to $s_i$ than to any other seed.  The boundary between two
adjacent cells is a straight **ridge**.  For a finite point set in the plane,
some ridges extend to infinity (open ridges); those are discarded.

### Why Tile Before Running Voronoi

`scipy.spatial.Voronoi` operates on an ordinary flat plane and has no concept
of periodic boundaries.  A seed point near the edge of the domain would have
an incomplete cell — the missing half is across the boundary.

The fix is to supply the algorithm with **ghost copies** of every seed: the
base cell is replicated into all 8 surrounding tiles (3×3 grid).  Extra seeds
now fill the neighbourhood of every boundary point and the Voronoi cells are
fully closed everywhere inside the base tile.  After the computation the ghost
ridges are discarded and only ridges with at least one seed in the base tile
are kept.

---

## Step-by-Step Walkthrough

### 1. Random Point Placement with Minimum Separation

```python
def _periodic_dist_sq(x1, y1, x2, y2, width, height):
    dx = min(abs(x1 - x2), width  - abs(x1 - x2))
    dy = min(abs(y1 - y2), height - abs(y1 - y2))
    return dx * dx + dy * dy
```

Distance on the torus in each axis is the smaller of the direct gap and the
wrap-around gap:

$$\delta x = \min\!\bigl(|x_1 - x_2|,\; W - |x_1 - x_2|\bigr)$$
$$\delta y = \min\!\bigl(|y_1 - y_2|,\; H - |y_1 - y_2|\bigr)$$
$$d^2 = \delta x^2 + \delta y^2$$

Returning the *squared* distance avoids a square-root and is sufficient for
comparisons against a threshold.

The generator loop uses rejection sampling:

```python
while len(points) < num_points and attempts < 1000:
    x, y = uniform(0, width), uniform(0, height)
    if all(_periodic_dist_sq(x, y, px, py, width, height) >= min_dist_sq
           for px, py in points):
        points.append((x, y))
        attempts = 0
    else:
        attempts += 1
```

A candidate point is accepted only if its periodic distance to every already-placed
point is at least $2r$ (the exclusion diameter).  The attempt counter resets on
each success; if 1000 consecutive candidates are rejected the domain is
considered full and the loop exits early.

The result is a **random hard-disk packing** on the torus: every point is at
least $2r$ from every other, respecting wrap-around.

---

### 2. Constructing the Tiled Point Set

```python
tile_offsets = [(dx * width, dy * height) for dx in range(-1, 2) for dy in range(-1, 2)]
```

This produces the 9 offsets:

| $(dx,dy)$ | Description         |
|-----------|---------------------|
| $(0, 0)$  | base cell           |
| $(\pm W, 0)$ | left / right copy |
| $(0, \pm H)$ | bottom / top copy |
| $(\pm W, \pm H)$ | corner copies |

Every base-cell point is replicated at each offset and labelled with its
base index and the tile shift:

```python
for i, (x, y) in enumerate(points):
    for dx, dy in tile_offsets:
        tiled_coords.append((x + dx, y + dy))
        tiled_meta.append((i, dx, dy))   # (base_index, shift_x, shift_y)
```

`tiled_meta[k]` makes it possible to look up, for any tiled point $k$, which
original seed it came from and by how much it was shifted.

---

### 3. Running the Voronoi Algorithm

```python
vor = Voronoi(np.array(tiled_coords))
```

`scipy.spatial.Voronoi` computes the complete Voronoi diagram for all
$9N$ tiled points.  The two key arrays used afterwards are:

| Array | Shape | Meaning |
|-------|-------|---------|
| `vor.ridge_points` | $(R, 2)$ | Indices of the two seed points on either side of each ridge |
| `vor.ridge_vertices` | $(R, 2)$ | Indices into `vor.vertices` for the two endpoints of the ridge; $-1$ = infinity |
| `vor.vertices` | $(V, 2)$ | Coordinates of all Voronoi vertices |

---

### 4. Edge Classification and Deduplication

The inner loop iterates over every ridge and classifies it:

```python
for (p, q), (vi, vj) in zip(vor.ridge_points, vor.ridge_vertices):
    if vi == -1 or vj == -1:   # open (infinite) ridge — discard
        continue
    bi, dxi, dyi = tiled_meta[p]
    bj, dxj, dyj = tiled_meta[q]
    is_base_p = dxi == 0 and dyi == 0
    is_base_q = dxj == 0 and dyj == 0
    if not (is_base_p or is_base_q):
        continue               # both seeds are ghosts — discard
```

**Open ridges** extend to infinity and carry no useful finite length.

**Both-ghost ridges** are pure artefacts of the tiling and are never needed.

The ridge length is computed only for ridges that pass the filter:

```python
length = math.dist(vor.vertices[vi], vor.vertices[vj])
```

#### Interior edge (both seeds in base tile)

```python
if is_base_p and is_base_q:
    key = (min(bi, bj), max(bi, bj))
    if key not in seen:
        edges.append((bi, bj, length))
        seen.add(key)
```

The canonical key `(min, max)` ensures the pair is stored once regardless of
which order the Voronoi algorithm emits it.

#### Periodic edge (one seed is a ghost)

```python
base, other, shift = (bi, bj, (dxj, dyj)) if is_base_p else (bj, bi, (dxi, dyi))
fwd = (base, other,  shift)
rev = (other, base, (-shift[0], -shift[1]))
if fwd not in seen and rev not in seen:
    periodic_edges.append((base, other, shift[0], shift[1], length))
    seen.add(fwd)
```

The ghost seed's tile offset is the **crossing shift**: when tiling the domain,
the neighbour `other` is reached from `base` by translating by `(shift_x, shift_y)`.

The forward / reverse deduplication is necessary because the same physical
boundary crossing appears *twice*:
- once with `base` in the base tile and `other`'s ghost at $+\text{shift}$,
- once (from the symmetric tile) with `other` in its own base tile and
  `base`'s ghost at $-\text{shift}$.

Checking both `fwd` and `rev` before inserting keeps exactly one record per
physical edge.

---

### 5. Adjacency Map

```python
adjacency = {i: [] for i in range(len(points))}
for i, j, _ in edges:
    adjacency[i].append((j,  0.0,  0.0))
    adjacency[j].append((i,  0.0,  0.0))
for i, j, dx, dy, _ in periodic_edges:
    adjacency[i].append((j,  dx,  dy))
    adjacency[j].append((i, -dx, -dy))
```

`adjacency[i]` is the list of neighbours of point $i$, each entry being
`(j, dx, dy)`.  The physical coordinates of neighbour $j$ as seen from $i$ are:

$$\bigl(x_j + dx,\; y_j + dy\bigr)$$

For interior edges $dx = dy = 0$.  For periodic edges the shift encodes how far
across the boundary the neighbour lives.

Both directions are populated: if point $i$ has neighbour $j$ with shift
$(dx, dy)$, then point $j$ has neighbour $i$ with shift $(-dx, -dy)$.  This
makes the adjacency map consistent for any algorithm that walks the graph
starting from an arbitrary node.

---

## Return Value

```python
return points, edges, periodic_edges, adjacency
```

| Name | Type | Content |
|------|------|---------|
| `points` | `list[(float, float)]` | Base-cell seed coordinates $(x, y)$ |
| `edges` | `list[(int, int, float)]` | Interior edges: `(i, j, length)` |
| `periodic_edges` | `list[(int, int, float, float, float)]` | Boundary-crossing edges: `(i, j, dx, dy, length)` |
| `adjacency` | `dict[int, list[(int, float, float)]]` | `i → [(j, dx, dy), …]` |

All indices refer to positions in `points`.

---

## Tiling the Result

To tile $n_x \times n_y$ copies of the base cell, a consumer would:

1. For each tile $(m, n)$ create points at $(x + m W,\; y + n H)$.
2. For interior edges use the duplicated indices directly.
3. For periodic edge $(i, j, dx, dy)$ in tile $(m, n)$, the partner point $j$
   lives in the tile at $(m + dx/W,\; n + dy/H)$ — the shift automatically
   gives the correct neighbour tile.

Because each periodic edge is stored once and the adjacency shift is
self-consistent, there are no duplicate edges and no missing connections at
tile boundaries.
