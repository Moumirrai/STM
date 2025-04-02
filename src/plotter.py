import matplotlib.pyplot as plt
from models import TrussData


def render_truss_structure(truss: TrussData):
    for node in truss.nodes.values():
        plt.plot(node.dx, node.dy, 'o', color='black')
        plt.text(node.dx, node.dy - 0.02, f"N{node.id}", fontsize=8, ha='right')

    for element in truss.elements:
        start_node = truss.nodes[element.starting_node]
        end_node = truss.nodes[element.ending_node]
        plt.plot([start_node.dx, end_node.dx], [start_node.dy, end_node.dy], 'k-')

    plt.xlabel('dx')
    plt.ylabel('dy')
    plt.title('Nodes and Elements')
    plt.grid(True)

    plt.gca().invert_yaxis()
    plt.show()