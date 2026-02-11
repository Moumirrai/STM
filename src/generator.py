from structure_parser import EigenstrainDefinition, ElementDefinition, NodeDefinition, StructureDefinition, MasterDefinition, DependencyDefinition


def create_cantilever_beam(length: float, height: float, nx: int, ny: int,
                          default_E: float = 210e6, default_A: float = 0.000004) -> StructureDefinition:
    """
    Create a simple cantilever beam structure.
    
    Fixed at left edge, free at right.
    """
    # Create nodes in a grid
    nodes = []
    dx_step = length / (nx - 1) if nx > 1 else 0
    dy_step = height / (ny - 1) if ny > 1 else 0
    
    for j in range(ny):
        for i in range(nx):
            dx = i * dx_step
            dy = j * dy_step
            
            # Fix left edge (i == 0)
            constraints = "xy" if i == 0 else ""
            
            node_def = NodeDefinition(
                dx=dx,
                dy=dy,
                constraints=constraints,
                deformations={},
                loads={}
            )
            nodes.append(node_def)
    
    # Create elements - all possible connections
    elements = []
    for j in range(ny):
        for i in range(nx):
            current = j * nx + i
            # Connect to right neighbor
            if i < nx - 1:
                right = j * nx + i + 1
                elements.append(ElementDefinition(
                    starting_node=current,
                    ending_node=right,
                    E=default_E,
                    A=default_A
                ))
            # Connect to bottom neighbor
            if j > 0:
                bottom = (j - 1) * nx + i
                elements.append(ElementDefinition(
                    starting_node=current,
                    ending_node=bottom,
                    E=default_E,
                    A=default_A
                ))
    
    volume = length * height
    
    return StructureDefinition(
        nodes=nodes,
        elements=elements,
        dependencies=[],
        eigenstrain=EigenstrainDefinition(x=0.0, y=0.0, angle=0.0),
        defaultYoungsModulus=default_E,
        defaultCrossSectionArea=default_A,
        volume=volume
    )
    
def create_tie_structure(x: float) -> StructureDefinition:
    width = 0.2
    height = 0.1
    
    x = height * x
    
    nodes = [
        NodeDefinition(dx=0, dy=0),
        NodeDefinition(dx=0, dy=height),
        NodeDefinition(dx=width, dy=0),
        NodeDefinition(dx=width, dy=height),
        
        NodeDefinition(dx=width/2, dy=height/2 - x/2),
        NodeDefinition(dx=width/2, dy=height/2 + x/2, constraints="xy"),
        
        NodeDefinition(dx=width/2, dy=height + height / 2 - x/2),
    ]
    
    elements = [
        ElementDefinition(starting_node=0, ending_node=1),
        ElementDefinition(starting_node=2, ending_node=3),
        ElementDefinition(starting_node=0, ending_node=4),
        ElementDefinition(starting_node=2, ending_node=4),
        ElementDefinition(starting_node=4, ending_node=5),
        ElementDefinition(starting_node=1, ending_node=5),
        ElementDefinition(starting_node=3, ending_node=5),
        ElementDefinition(starting_node=5, ending_node=6)
    ]
    
    dependencies = [
        DependencyDefinition(
            node=3,
            masters=[
                MasterDefinition(node=0, direction="x", factor=1.0),
                MasterDefinition(node=0, direction="y", factor=1.0)
            ]
        ),
        DependencyDefinition(
            node=2,
            masters=[
                MasterDefinition(node=0, direction="x", factor=1.0),
                MasterDefinition(node=0, direction="y", factor=1.0)
            ]
        ),
        DependencyDefinition(
            node=1,
            masters=[
                MasterDefinition(node=0, direction="x", factor=1.0),
                MasterDefinition(node=0, direction="y", factor=1.0)
            ]
        ),
        DependencyDefinition(
            node=6,
            masters=[
                MasterDefinition(node=4, direction="x", factor=1.0),
                MasterDefinition(node=4, direction="y", factor=1.0)
            ]
        )
    ]
    
    return StructureDefinition(
        nodes=nodes,
        elements=elements,
        dependencies=dependencies,
        eigenstrain=EigenstrainDefinition(x=1.0, y=0.0, angle=1.0),
    )
    