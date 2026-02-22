import matplotlib
import numpy as np
import math

from generator import create_tie_structure, create_tie_structure_angle
from parameter_solver import solveParameters_iso, solveParameters_orto

matplotlib.use("QtAgg")
import matplotlib.pyplot as plt

height = 0.1
width = height * 2
max_angle = math.degrees(math.atan(height / width))

x = np.linspace(0.001, max_angle, 100, endpoint=False)

Ex, Ey, vxy, vyx, Gxy = zip(*[solveParameters_orto(create_tie_structure_angle(height, width, i)) for i in x])

np.savetxt(
    "output.csv", np.column_stack((x, Ex, Ey, vxy, vyx, Gxy)), delimiter=",", comments=""
)

plt.plot(x, vxy)
plt.plot(x, vyx)
plt.show()
