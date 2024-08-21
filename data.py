#!/usr/bin/env python
# coding: utf-8
"""
Create Pandas DataFrame for C102
"""

import pandas as pd

# file path
file_path = "c102.txt"
# column name
columns = ["CUST NO.", "XCOORD.", "YCOORD.", "DEMAND", "READY TIME", "DUE DATE", "SERVICE TIME"]

def getData(file_path=file_path):
    # init data
    data = []
    # open file
    with open('c102.txt', 'r') as file:
        # skip first row
        next(file)
        # read file
        for line in file:
            if line.strip():
               values = line.split()
               data.append([int(values[0])] + [float(v) for v in values[1:]])
    # to DataFrame
    df = pd.DataFrame(data, columns=columns)
    return df

if __name__ == "__main__":
    df = getData()
    print(df.head())
