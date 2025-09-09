import json
import math
import random
from scipy.spatial import Delaunay

def generate_triangulation_json(x_points=6, y_points=5, x_range=(0, 8), y_range=(0, 10), fill_factor=0.5, border_margin=0.5):
    # Generate random 2D points
    points = []
    points.append([x_range[0], y_range[0]])
    points.append([x_range[0], y_range[1]])
    points.append([x_range[1], y_range[0]])
    points.append([x_range[1], y_range[1]])
    for _ in range(x_points-2):
        x = random.uniform(*x_range)
        points.append([x, y_range[0]])
        points.append([x, y_range[1]])
    for _ in range(y_points-2):
        y = random.uniform(*y_range)
        points.append([x_range[0], y])
        points.append([x_range[1], y])
    for _ in range(math.floor(x_range[1] * y_range[1] * fill_factor)):  # Add some random internal points
        x = random.uniform(x_range[0] + border_margin, x_range[1] - border_margin)
        y = random.uniform(y_range[0] + border_margin, y_range[1] - border_margin)
        points.append([x, y])

    print(points)
    
    # Perform Delaunay triangulation
    tri = Delaunay(points)
    
    # Nodes: list of points with dx, dy as strings
    nodes = [{"dx": str(p[0]), "dy": str(p[1])} for p in points]
    
    # Elements: edges from triangulation (unique edges)
    edges = set()
    for simplex in tri.simplices:
        for i in range(3):
            edge = tuple(sorted([int(simplex[i]), int(simplex[(i+1)%3])]))
            edges.add(edge)
    
    elements = [{"starting_node": int(e[0]), "ending_node": int(e[1])} for e in edges]
    
    # No dependencies or constraints for simplicity
    dependencies = []
    
    # Structure the JSON
    data = {
        "$schema": "./input.schema.json",
        "nodes": nodes,
        "dependencies": dependencies,
        "elements": elements
    }

    #save json to file
    with open("data/mesh.json", "w") as f:
        f.write(json.dumps(data, indent=2))

    return json.dumps(data, indent=2)

# Example usage
if __name__ == "__main__":
    generate_triangulation_json()