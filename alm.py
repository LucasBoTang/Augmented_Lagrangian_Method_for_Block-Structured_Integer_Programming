#!/usr/bin/env python
# coding: utf-8
"""
Augmented Lagrangian Method
"""

import numpy as np

import greedy
import utlis

def solve(df, num_customers=25, num_vehicles=3, k_max=100):
    """
    ALM main function
    """
    # get n customers
    df = df.iloc[:num_customers+1]
    # initialize problem
    x, λ, ρ = initialize(df, num_customers, num_vehicles)
    # get global constraints
    Aj, A, b = getGlobalConstraintsCoefficients(num_customers, num_vehicles)
    # iterations
    for k in range(k_max):
        # step size
        alpha = 0.1
        # update x
        x = updatePrimalSolution(df, num_customers, num_vehicles, x, λ, ρ, Aj, method="c")
        # constraints violation
        violation = A @ x.flatten() - b
        # termination condition
        if np.linalg.norm(violation) < 1e-3:
            return x, λ, ρ
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


def getGlobalConstraintsCoefficients(num_customers, num_vehicles):
    """
    get constraints coefficient A & b
    """
    Aj = utlis.generateAj(num_customers)
    A = np.hstack([Aj] * num_vehicles)
    b = np.ones(num_customers)
    return Aj, A, b


def updatePrimalSolution(df, num_customers, num_vehicles, x, λ, ρ, Aj, method):
    """
    update primal solution with BCD method
    """
    pass
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
