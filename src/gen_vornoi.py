import numpy as np
import re
from pathlib import Path
from vornoi_structure import generate_voronoi_structure
from parameter_solver import solveParameters_iso, solveParameters_orto
from models import TrussData
from plotter import plot_deformed_structure
from solver import TrussSolver
from structure_parser import parse_structure_data, truss_data_to_json_string


np.set_printoptions(
    linewidth=250,
)

trussData = generate_voronoi_structure(5, 5, 10000, 0.2, 2, 0.1)

truss: TrussData = parse_structure_data(
    trussData, explicitEigenStrain=np.array([0.0, 0.0, 0.0])
)

jsonString = truss_data_to_json_string(truss)

def get_next_structure_file_path() -> Path:
    structures_dir = Path(__file__).resolve().parent.parent / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    pattern = re.compile(r"^v(\d+)\.json$")
    indices = []
    for path in structures_dir.glob("v*.json"):
        match = pattern.match(path.name)
        if match:
            indices.append(int(match.group(1)))

    next_index = max(indices) + 1 if indices else 0
    return structures_dir / f"v{next_index}.json"


output_path = get_next_structure_file_path()
with output_path.open("w", encoding="utf-8") as f:
    f.write(jsonString)

print(f"Saved structure JSON to {output_path}")