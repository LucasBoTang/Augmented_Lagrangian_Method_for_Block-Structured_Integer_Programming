#!/usr/bin/env python
# coding: utf-8
"""
Utlity Functions
"""

import numpy as np

def generateAj(num_customers):
    # number of edges
    num_edges = num_customers * (num_customers + 1)
    # init Aj
    Aj = np.zeros((num_customers, num_edges))
    # fill the matrix
    for row in range(num_customers):
        s = row + 1
        Aj[row, s*num_customers:(s+1)*num_customers] =  1
    return Aj


if __name__ == "__main__":

    # generate Aj
    Aj = generateAj(3)
    print(Aj)
