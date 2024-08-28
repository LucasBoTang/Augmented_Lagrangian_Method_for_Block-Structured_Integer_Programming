# Augmented Lagrangian Method for Block-Structured Integer Programming

### Summary

The paper titled "A Customized Augmented Lagrangian Method for Block-Structured Integer Programming" investigates the problem of Block-Structured Integer Programming (BSIP) and proposes a customized Augmented Lagrangian Method (ALM) to solve it. They propose a new decomposition technique that breaks down the minimization of the augmented Lagrangian function into several simpler subproblems. The Block Coordinate Descent (BCD) method can solve these subproblems independently. By iteratively solving these subproblems and updating the Lagrange multipliers, the method converges towards a globally feasible and optimal solution.

### Code

We have reproduced the application of the Augmented Lagrangian Method to the Vehicle Routing Problem with Time Windows (CVRPTW).

In the repository, you can find the following code:

- `data.py`: Reads and processes the Solomon dataset.
- `alm.py`: Implements the core logic of the Augmented Lagrangian Algorithm.
- `calm.py`: Incorporates feasible solution finding to `alm.py`.
- `bcd.py`: Implements the Block Coordinate Descent method for updating decision variables in the Augmented Lagrangian Algorithm, supporting Proximal Linear and Classical Update strategies.
- `greedy.py`: Implements a greedy algorithm to find an initial feasible solution using the Block Coordinate Descent method.
- `feasible_heuristic.py`: Implements heuristic algorithms for finding feasible solutions during iterations, currently only includes the Sweeping Technique.
- `gurobi.py`: Models and solves the problem using integer programming with Gurobi.

### Report

For a detailed report on the reproduction process, please refer to the following report: [Report on Reproduction](https://github.com/LucasBoTang/Augmented_Lagrangian_Method_for_Block-Structured_Integer_Programming/blob/main/report.pdf).

### Citation

For more detailed information on the method and its application, please refer to the original paper:

Wang, R., Zhang, C., Pu, S., Gao, J., & Wen, Z. (2024). A Customized Augmented Lagrangian Method for Block-Structured Integer Programming. IEEE Transactions on Pattern Analysis and Machine Intelligence.

BibTeX:

```
@article{wang2024customized,
  title={A Customized Augmented Lagrangian Method for Block-Structured Integer Programming},
  author={Wang, Rui and Zhang, Chuwen and Pu, Shanwen and Gao, Jianjun and Wen, Zaiwen},
  journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
  year={2024},
  publisher={IEEE}
}
```
