#!/usr/bin/env python
# coding: utf-8
"""
Utlity Functions
"""

from scipy.spatial.distance import cdist
import numpy as np


def getCoefficients(df, num_customers, num_vehicles):
    """
    get constraints coefficient c, A & b
    """
    cj = generateCj(df)
    Aj = generateAj(num_customers)
    A = np.hstack([Aj] * num_vehicles)
    b = np.ones(num_customers)
    return cj, Aj, A, b


def generateCj(df):
    """
    get block of objective coefficient from distance matrix
    """
    # calculate distance matrix
    coords = df[["XCOORD.", "YCOORD."]].values
    distance_matrix = np.array(cdist(coords, coords, metric="euclidean"))
    # create a mask to remove diagonal elements
    mask = np.ones_like(distance_matrix, dtype=bool)
    np.fill_diagonal(mask, 0)
    # flatten the remaining elements
    cj = distance_matrix[mask].flatten()
    return cj


def generateAj(num_customers):
    """
    get block of global constraint coefficients
    """
    # number of edges
    num_edges = num_customers * (num_customers + 1)
    # init Aj
    Aj = np.zeros((num_customers, num_edges))
    # fill the matrix
    for row in range(num_customers):
        s = row + 1
        Aj[row, s*num_customers:(s+1)*num_customers] =  1
    return Aj


def computeGradient(x, cj, Aj, λ, ρ, violation):
    """
    Compute the gradient for the block j in the ALM.
    """
    grad_j = cj + Aj.T @ λ + ρ * (Aj.T @ violation)
    return grad_j


if __name__ == "__main__":

    # generate Aj
    Aj = generateAj(3)
    print(Aj)
