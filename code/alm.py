#!/usr/bin/env python
# coding: utf-8
"""
Augmented Lagrangian Method
"""

import time
import numpy as np
from tqdm import tqdm

import greedy
import utlis
import bcd

def solve(df, num_customers=25, num_vehicles=3, k_max=100, t_max=50, tol=1e-2, x_update_method="c"):
    """
    ALM main function
    """
    # get n customers
    df = df.iloc[:num_customers+1]
    # initialize problem
    x, λ, ρ = initialize(df, num_customers, num_vehicles)
    # get global constraints
    cj, Aj, c, A, b = utlis.getCoefficients(df, num_customers, num_vehicles)
    # initialize solution pool
    solution_pool = [[] for _ in range(num_vehicles)]
    # init step size
    α_0 = 1.0
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
        tqdm.write(f"Iteration {k+1}/{k_max}: Objective val = {objval:.4f}, Violation Norm = {violation_norm:.4f}")
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
    # no projection for equality constraints
    # λ = np.maximum(0, λ)
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
    #x, λ, ρ = solve(df, num_customers=50, num_vehicles=5)
