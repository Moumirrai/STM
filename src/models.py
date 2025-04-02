from dataclasses import dataclass
from typing import Dict, List, Optional
import numpy as np

@dataclass
class Node:
    id: int
    dx: float
    dy: float
    constrained_x: bool = False
    constrained_y: bool = False
    deformation_x: float = 0.0
    deformation_y: float = 0.0
    load_x: float = 0.0
    load_y: float = 0.0

@dataclass
class Element:
    starting_node: int
    ending_node: int
    sin: Optional[float] = None
    cos: Optional[float] = None
    length: Optional[float] = None

@dataclass
class TrussData:
    nodes: Dict[int, Node]
    elements: List[Element]