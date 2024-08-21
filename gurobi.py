#!/usr/bin/env python
# coding: utf-8
"""
Solve with Gurobi
"""
from scipy.spatial.distance import cdist
import numpy as np
import gurobipy as grb
from gurobipy import GRB

def solve(df, num_customers=25, num_vehicles=3):
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
    model = grb.Model("CVRP")
    # decision variables
    x = model.addVars(vehicles, edges, vtype=GRB.BINARY, name="x")
    w = model.addVars(vehicles, nodes, vtype=GRB.CONTINUOUS, name="w")
    # obj: min total dist
    model.setObjective(grb.quicksum(distance_matrix[s, t] * x[j, s, t]
                                    for j in vehicles
                                    for s, t in edges),
                       GRB.MINIMIZE)
    # constrs
    # visit constraint
    model.addConstrs(grb.quicksum(x[j, s, t] for j in vehicles for t in nodes
                                  if t != s) == 1
                     for s in nodes if s != 0)
    # flow balance constraint
    model.addConstrs(grb.quicksum(x[j, s, t]
                                  for t in nodes if t != s) ==
                     grb.quicksum(x[j, t, s]
                                  for t in nodes if t != s)
                     for j in vehicles for s in nodes)
    # departure and return constraint
    model.addConstrs(grb.quicksum(x[j, 0, t]
                                  for t in nodes if t != 0) == 1
                     for j in vehicles)
    model.addConstrs(grb.quicksum(x[j, s, 0]
                                  for s in nodes if s != 0) == 1
                     for j in vehicles)
    # capacity constraint
    model.addConstrs(grb.quicksum(demand[s] * x[j, s, t]
                                  for s, t in edges if s != t) <= capacity
                     for j in vehicles)
    # time window prerequisite constraint
    M = max(due_date)
    model.addConstrs(w[j, s] + (service_time[s] + distance_matrix[s, t]) -
                     M * (1 - x[j, s, t]) <= w[j, t]
                     for j in vehicles for s, t in edges if s != 0)
    # time window bound constraint
    model.addConstrs(w[j, s] >= ready_time[s] for j in vehicles for s in nodes)
    model.addConstrs(w[j, s] <= due_date[s] for j in vehicles for s in nodes)
    # solve
    model.optimize()
    return model


if __name__ == "__main__":

    import data

    # get data
    df = data.getData()

    # get model
    model = solve(df)
