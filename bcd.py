#!/usr/bin/env python
# coding: utf-8
"""
Block Coordinate Descent
"""

import utlis

def descent(df, num_customers, num_vehicles, x, λ, ρ, cj, Aj, A, b, method):
    # block coordinates descent
    for j in range(num_vehicles):
        # compute gradient
        grad_j = utlis.computeGradient(x, cj, Aj, λ, ρ, A, b)
        # update
    return x
