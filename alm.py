#!/usr/bin/env python
# coding: utf-8
"""
Augmented Lagrangian Method
"""

import numpy as np

import greedy
import utlis
import bcd

def solve(df, num_customers=25, num_vehicles=3, k_max=100, t_max=50, tol=1e-2):
    """
    ALM main function
    """
    # get n customers
    df = df.iloc[:num_customers+1]
    # initialize problem
    x, λ, ρ = initialize(df, num_customers, num_vehicles)
    # get global constraints
    cj, Aj, A, b = utlis.getCoefficients(df, num_customers, num_vehicles)
    # init step size
    alpha_0 = 0.1
    # iterations
    for k in range(k_max):
        # update x
        x = updatePrimalSolution(df, num_customers, num_vehicles,
                                 x, λ, ρ, cj, Aj, A, b,
                                 t_max, tol, method="c")
        # constraints violation
        violation = A @ x.flatten() - b
        # termination condition
        if np.linalg.norm(violation) < tol:
            return x, λ, ρ
        # update step size
        alpha = alpha_0 / np.sqrt(k+1)
        # update λ & ρ
        λ = updateLagrangeMultipliers(λ, violation, alpha)
        ρ = updatePenaltyCoefficient(ρ, σ=1.1)
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
                         x, λ, ρ, cj, Aj, A, b, t_max, tol, method):
    """
    update primal solution with BCD method
    """
    # iterations
    for t in range(t_max):
        x_new = bcd.descent(df, num_customers, num_vehicles,
                            x, λ, ρ, cj, Aj, A, b, method)
        # termination condition
        if np.linalg.norm(x_new - x) < tol:
            return x_new
        # update
        x = x_new
    return x


def updateLagrangeMultipliers(λ, violation, alpha):
    """
    update Lagrange multipliers via projected subgradient
    """
    # update with subgradient
    λ += alpha * violation
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
    x, λ, ρ = solve(df, num_customers=50, num_vehicles=5)
