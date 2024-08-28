# Review and Reproduction: Augmented Lagrangian Method for Block-Structured Integer Programming

We have reproduced the application of the Augmented Lagrangian Method to the Vehicle Routing Problem with Time Windows (CVRPTW).

In the repository, you can find the following code:

- `c102.txt`: Contains the C102 instance from the Solomon dataset.
- `data.py`: Reads and processes the Solomon dataset.
- `alm.py`: Implements the core logic of the Augmented Lagrangian Algorithm.
- `calm.py`: Incorporates feasible solution finding to `alm.py`.
- `bcd.py`: Implements the Block Coordinate Descent method for updating decision variables in the Augmented Lagrangian Algorithm, supporting Proximal Linear and Classical Update strategies.
- `greedy.py`: Implements a greedy algorithm to find an initial feasible solution using the Block Coordinate Descent method.
- `feasible_heuristic.py`: Implements heuristic algorithms for finding feasible solutions during iterations, currently only includes the Sweeping Technique.
- `gurobi.py`: Models and solves the problem using integer programming with Gurobi.
