#!/usr/bin/env python
# coding: utf-8
"""
Customized Augmented Lagrangian Method
"""

import time
import numpy as np
from tqdm import tqdm

import greedy
import feasible_heuristic
import utlis
import bcd

def solve(df, num_customers=25, num_vehicles=3, k_max=100, t_max=50, tol=1e-2,
          x_update_method="c", feasible_solution_method="s"):
    """
    ALM main function which find feasible solution during iterations
    """
    # get n customers
    df = df.iloc[:num_customers+1]
    # initialize problem
    x, λ, ρ = initialize(df, num_customers, num_vehicles)
    # get global constraints
    cj, Aj, c, A, b = utlis.getCoefficients(df, num_customers, num_vehicles)
    # initialize solution pool
    solution_pool = [[] for _ in range(num_vehicles)]
    best_solution = None
    best_obj_val = float("inf")
    optimality_gap = float("inf")
    # init step size
    α_0 = 0.1
    # init timer
    tick = time.time()
    # iterations
    for k in tqdm(range(k_max)):
        # update x
        x = updatePrimalSolution(df, num_customers, num_vehicles,
                                 x, λ, ρ, cj, Aj, A, b,
                                 t_max, tol, solution_pool, x_update_method)
        # constraints violation
        violation = A @ x.flatten() - b
        # violation norm
        violation_norm = np.linalg.norm(violation)
        # obj val
        objval = c @ x.flatten()
        # update tqdm description with violation norm
        tqdm.write(f"Iteration {k+1}/{k_max}: Objective val = {objval:.4f}, Violation Norm = {violation_norm:.4f}, Optimality Gap = {optimality_gap*100:.2f}%")
        # Find feasible solution if violation is detected
        if feasible_solution_method == "s":
            x_feasible = feasible_heuristic.sweeping(x, num_customers, num_vehicles, solution_pool, A, b)
        elif feasible_solution_method == "p":
            x_feasible = feasible_heuristic.packing(x, num_customers, num_vehicles, solution_pool)
        # update best feasible solution
        if x_feasible is not None:
            feasible_objval = c @ x_feasible.flatten()
            if feasible_objval < best_obj_val:
                best_obj_val = feasible_objval
                best_solution = x_feasible
        # calculate dual gap
        if best_solution is not None:
            optimality_gap = (best_obj_val - objval) / best_obj_val
        # termination condition
        if violation_norm < tol:
            tock = time.time()
            elpased = tock - tick
            time.sleep(1)
            print()
            print(f"Converged after {k+1} iterations within {elpased:.2f} sec")
            return x, λ, ρ
        # update step size
        α = α_0 / np.sqrt(k+1)
        # update λ & ρ
        λ = updateLagrangeMultipliers(λ, violation, α)
        ρ = updatePenaltyCoefficient(ρ, σ=1.1)
    # record time
    tock = time.time()
    elapsed = tock - tick
    time.sleep(1)
    print()
    print(f"Reached maximum iterations ({k_max}) after {elapsed:.2f} sec.")
    return x, λ, ρ


def initialize(df, num_customers, num_vehicles):
    """
    initialize ALM variables
    """
    x_init = greedy.solve(df, num_customers, num_vehicles)
    λ_init = np.zeros(num_customers)
    ρ_init = 1
    return x_init, λ_init, ρ_init


def updatePrimalSolution(df, num_customers, num_vehicles,
                         x, λ, ρ, cj, Aj, A, b, t_max, tol, solution_pool, method):
    """
    update primal solution with BCD method
    """
    # iterations
    for t in range(t_max):
        x_new = bcd.descent(df, num_customers, num_vehicles,
                            x, λ, ρ, cj, Aj, A, b, solution_pool, method)
        # termination condition
        if np.linalg.norm(x_new - x) < tol:
            return x_new
        # update
        x = x_new
    return x


def updateLagrangeMultipliers(λ, violation, α):
    """
    update Lagrange multipliers via projected subgradient
    """
    # update with subgradient
    λ += α * violation
    # project onto the non-negative orthant
    λ = np.maximum(0, λ)
    return λ


def updatePenaltyCoefficient(ρ, σ):
    """
    increase penalty coefficient
    """
    ρ *= σ
    return ρ


if __name__ == "__main__":

    import data

    # get data
    df = data.getData()

    # get model
    x, λ, ρ = solve(df)
    #x, λ, ρ = solve(df, num_customers=100, num_vehicles=10)
