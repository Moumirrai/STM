import matplotlib
import numpy as np
import math
import time

from generator import create_tie_structure, create_tie_structure_angle
from parameter_solver import solveParameters_iso, solveParameters_orto

#matplotlib.use("QtAgg")
import matplotlib.pyplot as plt

height = 0.1
width = height * 1
max_angle = math.degrees(math.atan(height / width))

x = np.linspace(0.001, max_angle, 50, endpoint=False)

results = []
total = len(x)
start_time = time.perf_counter()

for idx, angle in enumerate(x, start=1):
    elapsed = time.perf_counter() - start_time
    progress_message = f"Solving {idx}/{total} | elapsed: {elapsed:.2f}s"
    print(f"\r{progress_message:<80}", end="", flush=True)
    results.append(solveParameters_orto(create_tie_structure_angle(height, width, angle)))

total_elapsed = time.perf_counter() - start_time
final_message = f"Solved {total}/{total} | total: {total_elapsed:.2f}s"
print(f"\r{final_message:<80}")
Ex, Ey, vxy, vyx, Gxy = zip(*results)

np.savetxt(
    "output.csv", np.column_stack((x, Ex, Ey, vxy, vyx, Gxy)), delimiter=",", comments=""
)

plt.plot(x, vxy, label="vxy")
plt.plot(x, vyx, label="vyx")
plt.xlabel("Angle")
plt.legend()
plt.show()
