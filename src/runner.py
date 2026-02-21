import matplotlib
import numpy as np

from generator import create_tie_structure
from parameter_solver import solveParameters_iso, solveParameters_orto

matplotlib.use("QtAgg")
import matplotlib.pyplot as plt

x = np.linspace(0.001, 1, 100, endpoint=False)

Ex, Ey, vxy, vyx, Gxy = zip(*[solveParameters_orto(create_tie_structure(i)) for i in x])

np.savetxt(
    "output.csv", np.column_stack((Ex, Ey, vxy, vyx, Gxy)), delimiter=",", comments=""
)

plt.plot(x, vxy)
plt.plot(x, vyx)
plt.show()
