#!/usr/bin/env python
# coding: utf-8
"""
Hueristic Algorithm to find feasible solution
"""

from copy import deepcopy
import numpy as np

def sweeping(x, num_customers, num_vehicles, solution_pool, A, b, q_max=10, tol=1e-2):
    """
    sweeping technique to find a feasible solution
    """
    # deepo copy
    solution_pool = deepcopy(solution_pool)
    # init feasible solution
    x_feasible = np.zeros_like(x)
    for j in range(num_vehicles):
        x_feasible[j] = solution_pool[j].pop()
    for _ in range(q_max):
        # iterate through all blockss
        for j in range(num_vehicles):
            #print(np.sum(x_feasible))
            # solution pool is empty
            if not solution_pool[j]:
                candidate = np.zeros_like(x[j])
            else:
                # get a candidate from solution pool
                candidate = solution_pool[j].pop()
            # check feasibility
            if checkFeasibility(x_feasible, j, candidate, A, b, tol):
                # get feasible solution
                x_feasible[j] = candidate
                return x_feasible
            else:
                # assign 0
                x_feasible[j] = np.zeros_like(x[j])
    # no feasible solution
    return None


def checkFeasibility(x, j, candidate, A, b, tol):
    """
    check feasibility of global constraints A x = b for current candidate
    """
    # deep copy
    x = deepcopy(x)
    # assign candidate
    x[j] = candidate
    # constraints violation
    violation = np.abs(A @ x.flatten() - b)
    return np.all(violation <= tol)


def packing(x, num_customers, num_vehicles, solution_pool):
    """
    packing technique to find a feasible solution
    """
    pass
