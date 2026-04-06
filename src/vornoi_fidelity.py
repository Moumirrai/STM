from concurrent.futures import ProcessPoolExecutor, as_completed
from time import perf_counter

import numpy as np

from parameter_solver import solveParameters_iso
from vornoi_structure import generate_voronoi_structure

num_tries = 1000
x = np.logspace(np.log10(0.6), np.log10(0.05), 25)  # logaritmicke rozlozeni

print(x)


def clear_terminal():
    # ANSI clear + cursor home; works on Linux terminals.
    print("\033[2J\033[H", end="")


def single_run(x_i):
    trussData = generate_voronoi_structure(5, 5, 50000, x_i)
    return solveParameters_iso(trussData)


def main():
    datapoints = []
    start_time = perf_counter()
    total_steps = len(x)

    for step_idx, x_i in enumerate(x, start=1):
        step_results = []

        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(single_run, x_i) for _ in range(num_tries)]
            for iteration_idx, future in enumerate(as_completed(futures), start=1):
                step_results.append(future.result())
                elapsed = perf_counter() - start_time

                clear_terminal()
                print(f"Step: {step_idx}/{total_steps}")
                print(f"Iteration: {iteration_idx}/{num_tries}")
                print(f"x: {x_i:.6f}")
                print(f"Elapsed from start: {elapsed:.2f} s")

        E, v = zip(*step_results)
        datapoints.append((x_i, np.mean(E), np.std(E), np.mean(v), np.std(v)))

    np.savetxt("vornoi_fidelity.csv", np.array(datapoints), delimiter=",", comments="")
    total_elapsed = perf_counter() - start_time
    clear_terminal()
    print(f"Done. Saved vornoi_fidelity.csv in {total_elapsed:.2f} s")


if __name__ == "__main__":
    main()