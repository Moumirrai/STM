import numpy as np
from scipy.optimize import least_squares
from termcolor import colored

from generator import create_cantilever_beam, create_tie_structure
from models import TrussData
from plotter import export_vtk
from solver import TrussSolver
from structure_parser import StructureDefinition, parse_json_file, parse_structure_data

np.set_printoptions(
    linewidth=250,
)


eigenstrainSets = [
    np.array([1, 0, 0]),
    np.array([0, 1, 0]),
    np.array([0, 0, 1]),
]


def solveParameters_iso(structure: StructureDefinition):
    results = []
    for eigenstrain in eigenstrainSets:
        truss: TrussData = parse_structure_data(
            structure, explicitEigenStrain=eigenstrain
        )
        solver = TrussSolver(truss)
        res = solver.solve()
        results.append(res)

    Ds = np.array(results).T

    def compute_D(params):
        E, v = params
        factor = E / (1 - v**2)
        return factor * np.array([[1, v, 0], [v, 1, 0], [0, 0, (1 - v) / 2]])

    def residuals(params):
        D_pred = compute_D(params)
        diff = (D_pred - Ds).flatten()
        # we divide the difference by the norm of Ds to normalize the residuals for better numerical stability and cost scale
        return diff / (
            np.linalg.norm(Ds) + 1e-10
        )  # add small value to avoid division by zero
        # return diff # same results without cost normalization

    # initial guesses and bounds
    initial_guess = [210e6, 0.4]  # E in Pa, v dimensionless
    bounds = ([1e2, 0], [1e12, 0.5])

    result = least_squares(residuals, initial_guess, bounds=bounds)
    fitted_E, fitted_v = result.x

    print(colored(f"D matrix from DOF elimination solver:\n{Ds}\n", "cyan"))

    D_fitted = compute_D([fitted_E, fitted_v])
    print(colored(f"Fitted D matrix:\n{D_fitted}\n", "light_green"))

    print(
        colored(
            f" E = {fitted_E:.2e} Pa \t v = {fitted_v:.3f} ", "black", "on_light_green"
        )
    )
    print("")
    # print("Optimization success:", result.success)
    print(colored(f"Final cost: {result.cost}", "light_red"))

    return fitted_E, fitted_v


def solveParameters_orto(structure: StructureDefinition):
    results = []
    for eigenstrain in eigenstrainSets:
        truss: TrussData = parse_structure_data(
            structure, explicitEigenStrain=eigenstrain
        )
        solver = TrussSolver(truss)
        res = solver.solve()
        results.append(res)

    Ds = np.array(results).T

    print(colored(f"D matrix from DOF elimination solver:\n{Ds}\n", "cyan"))

    def compute_D(params):
        Ex, Ey, vxy, vyx, Gxy = params
        factor = 1 / (1 - vxy * vyx)
        return factor * np.array(
            [
                [Ex, vyx * Ex, 0],
                [vxy * Ey, Ey, 0],
                [0, 0, Gxy * (1 - vxy * vyx)],
            ]  # division may be off
        )

    def residuals(params):
        D_pred = compute_D(params)
        diff = (D_pred - Ds).flatten()
        # we divide the difference by the norm of Ds to normalize the residuals for better numerical stability and cost scale
        return diff / (
            np.linalg.norm(Ds) + 1e-10
        )  # add small value to avoid division by zero
        # return diff # same results without cost normalization

    # initial guesses and bounds
    initial_guess = [210e6, 210e6, 0.4, 0.4, 21e6]  # E in Pa, v dimensionless
    bounds = ([1e2, 1e2, 0, 0, 0], [1e12, 1e12, 0.5, 0.5, 1e12])

    result = least_squares(residuals, initial_guess, bounds=bounds)
    f_Ex, f_Ey, f_vxy, f_vyx, f_Gxy = result.x

    D_fitted = compute_D([f_Ex, f_Ey, f_vxy, f_vyx, f_Gxy])

    return f_Ex, f_Ey, f_vxy, f_vyx, f_Gxy
