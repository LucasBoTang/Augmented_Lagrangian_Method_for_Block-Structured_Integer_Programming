#!/usr/bin/env python
# coding: utf-8
"""
Augmented Lagrangian Method
"""

import numpy as np

import greedy

def solve(df, num_customers=25, num_vehicles=3):
    """
    ALM main function
    """
    # get n customers
    df = df.iloc[:num_customers+1]
    # initialize problem
    x, λ, ρ = initialize(df, num_customers, num_vehicles)


def initialize(df, num_customers, num_vehicles):
    """
    initialize ALM variables
    """
    x_init = greedy.solve(df, num_customers, num_vehicles)
    λ_init = np.zeros(num_vehicles)
    ρ_init = 1
    return x_init, λ_init, ρ_init


if __name__ == "__main__":

    import data

    # get data
    df = data.getData()

    # get model
    model = solve(df)
