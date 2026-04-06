import matplotlib
import numpy as np
import math
import time
from concurrent.futures import ProcessPoolExecutor

from vornoi_structure import generate_voronoi_structure, aniotropise
from parameter_solver import solveParameters_iso, solveParameters_orto

#matplotlib.use("QtAgg")
import matplotlib.pyplot as plt

def _run_single_try(iter_val):
    """Helper function for multiprocessing: generates structure and solves for one try."""
    return solveParameters_orto(aniotropise(generate_voronoi_structure(6, 6, 10000, 0.2), 100e9, iter_val))

if __name__ == "__main__":
    x = np.linspace(1e9, 100e9, 70, endpoint=True)

    results = []
    total = len(x)
    start_time = time.perf_counter()

    for idx, iter_val in enumerate(x, start=1):
        elapsed = time.perf_counter() - start_time
        progress_message = f"Solving {idx}/{total} | elapsed: {elapsed:.2f}s"
        print(f"\r{progress_message:<80}", end="", flush=True)
        
        # Run 30 tries in parallel across all cores
        step_results = []
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(_run_single_try, iter_val) for _ in range(500)]
            step_results = [f.result() for f in futures]
        
        Ex_values, Ey_values, vxy_values, vyx_values, Gxy_values = zip(*step_results)
        avg_Ex = np.mean(Ex_values)
        avg_Ey = np.mean(Ey_values)
        avg_vxy = np.mean(vxy_values)
        avg_vyx = np.mean(vyx_values)
        avg_Gxy = np.mean(Gxy_values)
        results.append((avg_Ex, avg_Ey, avg_vxy, avg_vyx, avg_Gxy))

    total_elapsed = time.perf_counter() - start_time
    final_message = f"Solved {total}/{total} | total: {total_elapsed:.2f}s"
    print(f"\r{final_message:<80}")
    Ex, Ey, vxy, vyx, Gxy = zip(*results)

    np.savetxt(
        "vornoi_out.csv", np.column_stack((x, Ex, Ey, vxy, vyx, Gxy)), delimiter=",", comments=""
    )

    plt.plot(x, Ex, label="Ex")
    plt.plot(x, Ey, label="Ey")
    plt.xlabel("Angle")
    plt.legend()
    plt.show()
