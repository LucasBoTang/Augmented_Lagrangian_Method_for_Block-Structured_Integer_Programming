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

### Model

The mathematical model for the Capacitated Vehicle Routing Problem with Time Windows (CVRPTW) is formulated as follows:

$$
\begin{align}
\min & \sum_{j \in N_p} \sum_{(s,t) \in E} d_{st} x^j_{st} \\
\text{s.t.} \quad & \sum_{j \in N_p} \sum_{t \in V : t \neq s} x^j_{st} = 1, \quad \forall s \in V \setminus \{0\}, \\
& \sum_{t \in V \setminus \{s\}} x^j_{st} = \sum_{t \in V \setminus \{s\}} x^j_{ts}, \quad \forall s \in V, \forall j \in N_p, \\
& \sum_{t \in V \setminus \{0\}} x^j_{0t} = 1, \quad \forall j \in N_p, \\
& \sum_{s \in V} \sum_{t \in V \setminus \{s\}} c_s x^j_{st} \leq C, \quad \forall j \in N_p, \\
& w^j_s + T_{st} - M(1 - x^j_{st}) \leq w^j_t, \quad \forall (s, t) \in E, \forall j \in N_p, \\
& a_s \leq w^j_s \leq b_s, \quad \forall s \in V, \forall j \in N_p, \\
& x^j_{st} \in \{0, 1\}, \quad \forall (s, t) \in E, \forall j \in N_p, \\
& w^j_s \geq 0, \quad \forall s \in V, \forall j \in N_p.
\end{align}
$$

Let λ be the Lagrange multiplier associated with constraint (12), and ρ > 0 be the coefficient of the augmented term. The augmented Lagrangian function can then be expressed as:

$$
f(\mathbf{x}) + \sum_{s \in V \setminus \{0\}} \lambda_s \left( \sum_{j \in N_p} \sum_{t \in V : t \neq s} x^j_{st} - 1 \right) + \frac{\rho}{2} \sum_{s \in V \setminus \{0\}} \left( \sum_{j \in N_p} \sum_{t \in V : t \neq s} x^j_{st} - 1 \right)^2
$$

where $f(\mathbf{x})$ is the original objective function.

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

### License

This project is licensed under the MIT License - see the [LICENSE file](https://github.com/LucasBoTang/Augmented_Lagrangian_Method_for_Block-Structured_Integer_Programming/blob/main/LICENSE) for details.
