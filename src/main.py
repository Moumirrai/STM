from models import TrussData
from solver import TrussSolver
from plotter import export_vtk
from structure_parser import parse_json_file
from numpy import set_printoptions

set_printoptions(
    linewidth=250,
)


truss: TrussData = parse_json_file('./data/simple_truss1.json')

solver = TrussSolver(truss)

solver.solve()

export_vtk(truss)

exit(0)
