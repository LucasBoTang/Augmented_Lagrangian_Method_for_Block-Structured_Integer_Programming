#!/usr/bin/env python
# coding: utf-8
"""
Block Coordinate Descent
"""

from scipy.spatial.distance import cdist
import numpy as np
import gurobipy as grb
from gurobipy import GRB

import utlis

def descent(df, num_customers, num_vehicles, x, λ, ρ, cj, Aj, A, b, method):
    """
    Block coordinates descent to update xj
    """
    # block coordinates descent
    for j in range(num_vehicles):
        # update xj
        if method == "c":
            # classical update
            x[j] = classicalUpdate(df, num_customers, num_vehicles, x, j, cj, Aj, λ, ρ)
        if method == "p":
            x[j] = proximalLinearUpdate(df, num_customers, x, j, cj, Aj, A, b, λ, ρ)
    return x


def classicalUpdate(df, num_customers, num_vehicles, x, j, cj, Aj, λ, ρ):
    """
    classical update step for block j using Gurobi to solve the integer linear problem
    """
    # number of nodes
    num_nodes = num_customers + 1
    # get sets
    edges = [(i, j) for i in range(num_nodes) for j in range(num_nodes) if i != j]
    vehicles = list(range(num_vehicles))
    nodes = list(range(num_nodes))
    # get n customers
    df = df.iloc[:num_nodes]
    # calculate distance matrix
    coords = df[["XCOORD.", "YCOORD."]].values
    distance_matrix = np.array(cdist(coords, coords, metric="euclidean"))
    # get values for modeling
    demand = df['DEMAND'].values
    ready_time = df['READY TIME'].values
    due_date = df['DUE DATE'].values
    service_time = df['SERVICE TIME'].values
    # set capacity
    capacity = 200
    # init model
    model = grb.Model("classical_update")
    # disable output
    model.setParam('OutputFlag', 0)
    # decision variables
    xj = model.addVars(edges, vtype=GRB.BINARY, name="x")
    wj = model.addVars(nodes, vtype=GRB.CONTINUOUS, name="w")
    # obj
    obj_coeffs = cj + Aj.T @ λ
    temp = np.zeros_like(λ)
    for l in vehicles:
        if l != j:
            temp += Aj @ x[l]
    obj_coeffs += ρ * Aj.T @ (temp - 1/2)
    model.setObjective(grb.quicksum(obj_coeffs[i] * xj[s, t]
                                    for i, (s, t) in enumerate(edges)),
                       GRB.MINIMIZE)
    # constrs
    # flow balance constraint
    model.addConstrs(grb.quicksum(xj[s, t]
                                  for t in nodes if t != s) ==
                     grb.quicksum(xj[t, s]
                                  for t in nodes if t != s)
                     for s in nodes)
    # departure and return constraint
    model.addConstr(grb.quicksum(xj[0, t]
                                  for t in nodes if t != 0) == 1)
    model.addConstr(grb.quicksum(xj[s, 0]
                                  for s in nodes if s != 0) == 1)
    # capacity constraint
    model.addConstr(grb.quicksum(demand[s] * xj[s, t]
                                  for s, t in edges if s != t) <= capacity)
    # time window prerequisite constraint
    M = max(due_date)
    model.addConstrs(wj[s] + (service_time[s] + distance_matrix[s, t]) -
                     M * (1 - xj[s, t]) <= wj[t]
                     for s, t in edges if s != 0)
    # time window bound constraint
    model.addConstrs(wj[s] >= ready_time[s] for s in nodes)
    model.addConstrs(wj[s] <= due_date[s] for s in nodes)
    # solve
    model.optimize()
    # to numpy
    xj = utlis.sol2Numpy(xj)
    return xj


def proximalLinearUpdate(df, num_customers, x, j, cj, Aj, A, b, λ, ρ):
    """
    """
    pass
    return x[j]


def toNumpy(xj):
    """
    convert
    """
    # init
    xj_array = np.zeros(len(xj))
    # assign value
    for i, e in enumerate(xj):
        xj_array[i] = xj[e].X
    return xj_array
