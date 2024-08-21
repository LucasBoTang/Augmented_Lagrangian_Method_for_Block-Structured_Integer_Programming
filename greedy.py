#!/usr/bin/env python
# coding: utf-8
"""
Greedy Algorithm for CVRP
"""

from scipy.spatial.distance import cdist
import numpy as np

def solve(df, num_customers=25, num_vehicles=3):
    """
    Greedy algorithm which is block feasible
    """
    # get n customers
    df = df.iloc[:num_customers+1]
    # init
    vehicle_routes = []
    remaining_customers = set(range(1,num_customers+1))
    depot = df.iloc[0]
    vehicles_used = 0
    vehicle_capacity = 200
    # calculate distance matrix
    coords = df[["XCOORD.", "YCOORD."]].values
    distance_matrix = np.array(cdist(coords, coords, metric="euclidean"))
    while remaining_customers and vehicles_used < num_vehicles:
        # start from depot
        current_route = [0]
        current_load = 0
        current_time = depot["READY TIME"]
        current_position = 0
        while remaining_customers:
            # distance from current current_position
            distances = distance_matrix[current_position]
            # find earlest due date
            sorted_customers = sorted(remaining_customers, key=lambda x: df.loc[x]["DUE DATE"] - distance_matrix[current_position, x])
            found_customer = False
            # check validation
            for customer in sorted_customers:
                demand = df.loc[customer]["DEMAND"]
                ready_time = df.loc[customer]["READY TIME"]
                due_date = df.loc[customer]["DUE DATE"]
                service_time = df.loc[customer]["SERVICE TIME"]
                if current_load + demand <= vehicle_capacity and current_time + distances[list(remaining_customers).index(customer)] <= due_date:
                    current_route.append(customer)
                    current_load += demand
                    current_time += distances[list(remaining_customers).index(customer)] + service_time
                    current_position = customer
                    remaining_customers.remove(customer)
                    found_customer = True
                    break
            # no available customer
            if not found_customer:
                break
        # return to depot
        current_route.append(0)
        vehicle_routes.append(current_route)
        vehicles_used += 1
    # transform solution
    x = transformSolution(vehicle_routes, num_customers)
    return x


def transformSolution(vehicle_routes, num_customers):
    """
    transform from vehicle_routes to solution vector
    """
    # init
    edges = [(i, j) for i in range(num_customers+1) for j in range(num_customers+1) if i != j]
    x = np.zeros((len(vehicle_routes), len(edges)))
    # dict: edge to index
    e2i_table = {}
    for i, e in enumerate(edges):
        e2i_table[e] = i
    # route for each vehicle
    for j, route in enumerate(vehicle_routes):
        for s, t in zip(route[:-1], route[1:]):
            i = e2i_table[(s,t)]
            x[j,i] = 1
    return x

if __name__ == "__main__":

    import data

    # get data
    df = data.getData()

    # get model
    x = solve(df)
