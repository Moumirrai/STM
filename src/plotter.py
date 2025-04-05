import matplotlib.pyplot as plt
from models import TrussData


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