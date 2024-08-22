# 论文评论及复现：分块结构整数规划的增广拉格朗日方法

### 1. 概要

“A Customized Augmented Lagrangian Method for Block-Structured Integer Programming”这篇文章研究了具有分块结构整数线性规划（Block-Structured Integer Programming, BSIP）问题，并提出了一种定制化的增广拉格朗日方法来求解这类问题。文章中提到，这类问题因其在列车时刻表、车辆路径规划等众多实际应用中的普遍性而备受关注。但由于这类问题存在整数变量，是NP-hard的，不存在已知的多项式时间算法来提供精确解。

因此，作者在传统的增广拉格朗日法(Augmented Lagrangian Method, ALM)迭代求解的基础上进行了创新，通过引入一种新的分解技术，将增广拉格朗日函数的最小化分解为多个简单的子问题。这些子问题可以独立地使用块坐标下降法（Block Coordinate Descent, BCD）来求解。作者也在理论上证明了始问题和增广拉格朗日对偶问题之间的强对偶性，也对块坐标下降方法和增广拉格朗日方法的收敛性进行了理论分析。

### 2. 分块结构整数规划问题模型 

我们考虑一个带有块结构的整数规划问题，该问题的数学模型可以表示为：

$$ 
\begin{align}
f^{\text{IP}} := \min_{\mathbf{x}} \quad & \mathbf{c}^{\intercal} \mathbf{x} \\
\text{s.t.} \quad & \mathbf{A} \mathbf{x} \leq \mathbf{b}, \\
& \mathbf{x}_j \in \mathcal{X}_j, \quad j = 1, 2, \ldots, p,
\end{align}
$$

其中：

- $\mathbf{x}$是决策变量向量，它被分为$p$个子块。
- $\mathbf{c}^{\intercal} \mathbf{x}$是线性目标函数，$\mathbf{c}$是目标函数的系数向量。
- $\mathbf{A} \mathbf{x} \leq \mathbf{b}$是全局线性不等式约束，这些约束将所有的决策变量 $\mathbf{x}$ 耦合在一起。
- 每个子块$\mathbf{x}_j$需要满足特定的整数约束条件$\mathcal{X}_j = \{\mathbf{x}_j \in \{0,1\}^{n_j}: \mathbf{B}_j \mathbf{x}_j \leq \mathbf{D}_j\}$。

### 3. 增广拉格朗日方法

拉格朗日增广是一种用于解决带约束的优化问题的策略，它在拉格朗日函数的基础上，增加约束的二次惩罚项。以下是一个分块结构整数规划的拉格朗日增广模型：
$$L(\mathbf{x}, \mathbf{\lambda}, \rho) = \mathbf{c}^{\intercal} \mathbf{x} + \mathbf{\lambda}^{\intercal} (\mathbf{A} \mathbf{x} - \mathbf{b}) + \frac{\rho}{2} \| (\mathbf{A} \mathbf{x} - \mathbf{b})_{+} \|^2$$

其中，$\mathbf{\lambda} \geq \mathbf{0}$是拉格朗日乘子向量，$\rho > 0$是二次项的惩罚参数。相应的增广拉格朗日松弛问题如下：
$$
d(\mathbf{\lambda}, \rho) := \min_{\mathbf{x} \in \mathcal{X}} L(\mathbf{x}, \mathbf{\lambda}, \rho)
$$

我们称以下最大化问题为增广拉格朗日对偶问题:
$$
f_{\rho}^{\text{LD}} :=  \max_{\mathbf{\lambda} > \mathbf{0}} d(\mathbf{\lambda}, \rho)
$$
该对偶问题求解的是增广拉格朗日函数的最小值，并通过调整拉格朗日乘子$\mathbf{\lambda}$来最大化这个最小值。

#### 3.1 强对偶性

作者证明，如果原问题可行且其最优值有界，存在一个有限的惩罚参数值$\rho^*$，当$\forall \rho \geq \rho^*$，原问题$f^{\text{IP}}$和增广拉格朗日对偶问题$f_{\rho}^{\text{LD}}$之间存在强对偶性——它们的最优值相等。因此，我们可以通过解决对偶问题来找到原始问题的最优解。

#### 3.2 迭代步骤

增广拉格朗日方法的迭代过程通常包括以下步骤：

1. 在每一步 $k$，我们需要求解子问题以更新决策变量$\mathbf{x}^{k+1} := \arg \min_{\mathbf{x} \in \mathcal{X}} L(\mathbf{x}, \mathbf{\lambda}^k, {\rho}^k)$
2. 更新乘子$\mathbf{\lambda}^{k+1}$
3. 更新惩罚参数${\rho}^{k+1}$

由于增广拉格朗日松弛函数 $d(\mathbf{\lambda}, \rho)$ 是一个不可微的凹函数，我们需要利用当前迭代点 $\mathbf{x}^k$ 来构造该函数在 $(\mathbf{\lambda}^k, {\rho}^k)$ 处的次梯度，以便更新乘子 $\mathbf{\lambda}$ 和惩罚参数 $\rho$：
$$
\mathbf{A} \mathbf{x}^k - \mathbf{b} \in \partial_{\mathbf{\lambda}} d(\mathbf{\lambda}^k, {\rho}^k),
$$
$$
\frac{1}{2} \| (\mathbf{A} \mathbf{x}^k - \mathbf{b})_+ \| \in \partial_{\rho} d(\mathbf{\lambda}^k, {\rho}^k)
$$

鉴于乘子 $\mathbf{\lambda}$ 和惩罚参数 $\rho$ 必须满足非负条件，即$\mathbf{\lambda} \geq \mathbf{0}, \rho > 0$，我们可以通过投影次梯度方法（Projected Subgradient Method）来更新这些参数，确保它们在更新后仍然满足非负性约束：

$$
\begin{align}
\mathbf{\lambda}^{k+1} & := (\mathbf{\lambda}^k + \mathbf{\alpha}^k (\mathbf{A} \mathbf{x}^k - \mathbf{b}))_+, \\
{\rho}^{k+1} & := {\rho}^{k} + \frac{\mathbf{\alpha}^k}{2} \| (\mathbf{A} \mathbf{x}^k - \mathbf{b})_+ \|。
\end{align}
$$

然而，在论文的实验部分更新参数的方式有所不同。具体来说，作者使用了比例因子$\sigma$逐步增加惩罚参数$\rho$，以加大对违反约束条件的惩罚力度。因此，$\rho$实际上会按以下方式更新：
$$
\rho^{k+1} := \sigma  \cdot \rho^{k}
$$

#### 终止条件

算法在每次迭代后会检查当前解是否满足原始问题的约束条件，即：
$$
\| (\mathbf{A} \mathbf{x}^k - \mathbf{b})_+ \| = 0
$$
当这个值为零时，说明当前解最终的解在可行域内，且满足所有原始问题的约束，对偶问题已经收敛于$f_{\rho}^{\text{LD}} = f^{\text{IP}}$，因此可以终止算法。

#### 3.3 初始点的选择

尽管论文中并未详细说明，在增广拉格朗日方法中，初始点的选择将遵循以下规则：

- 初始决策变量$\mathbf{x}^0$需要在可行域$\mathcal{X}$内，以确保块坐标下降（BCD）方法能够执行并在有限次迭代后得到问题的块最优解。

- 初始拉格朗日乘子$\mathbf{\lambda}^0$应在非负域$\mathbb{R}_+^m$中。这些乘子通常初始化为零向量，即$\mathbf{\lambda}^0 = \mathbf{0}$。

- 初始惩罚参数$\rho$通常设定为一个较小的正值，并在迭代过程中逐渐增大，以保证收敛和满足约束条件。

#### 3.4 步长和收敛性

在增广拉格朗日方法中，步长 $\mathbf{\alpha}^k$ 控制着拉格朗日乘子的更新速率。通过选择合适的步长，可以使算法在经过有限次迭代后，使得对偶函数的值接近于原问题的最优值。

步长的计算将基于对偶函数$d(\mathbf{\lambda}^k, {\rho}^k)$次梯度的模长，按照如下公式确定：
$$
\mathbf{\alpha}^k = \frac{{\beta}^k}{\| \nabla d(\mathbf{\lambda}^k, {\rho}^k) \|}。
$$

作者通过理论分析，提出了${\beta}^k$的一个收敛条件：${\beta}^k = \sqrt{\frac{\theta}{K}}$，并给出了该收敛条件下对偶间隙的上界
$$\frac{\max_{\mathbf{x} \in \mathcal{x}} {\| \mathbf{A} \mathbf{x} - \mathbf{b} \|}^4}{2} \sqrt{\frac{5 \theta}{K}}。
$$
其中，$K$是最大迭代次数，$\theta$是初始解$(\mathbf{\lambda}^0, {\rho}^0)$距离最优解集合$S$的最小距离。

为了确保算法最终能够收敛到原问题的最优解 $(\mathbf{\lambda}^*, {\rho}^*)$，文中进一步提出了两个关键条件：

- 发散条件：确保算法能够无限次迭代直到收敛。$$\sum_{k=1}^{\infty} \beta_k = +\infty$$
- 平方和收敛条件：保证步长序列的平方和是有限的，有助于确保算法的稳定性。$$\sum_{k=1}^{\infty} \beta_k^2 < +\infty$$

尽管理论上提供了步长系数 ${\beta}^k$ 的选择依据，但文中并未具体说明在实际应用中应如何确定 ${\beta}^k$。在实践中，${\beta}^k$ 可能是一个正的递减序列，这样的设计旨在确保步长既不会过大以免引起剧烈更新，也不会过小以避免算法收敛过慢。


### 4. 块坐标下降法更新决策变量

在增广拉格朗日方法的迭代过程中，决策变量$\mathbf{x} := \arg \min_{\mathbf{x} \in \mathcal{X}} L(\mathbf{x}, \mathbf{\lambda}^k, {\rho}^k)$没有封闭形式的解，但我们可以应用迭代方法来精确或近似地求解它。

增广拉格朗日函数通过将全局约束嵌入目标函数，实现了约束条件的解耦。然而，由于函数中包含的二次项$\frac{\rho}{2} \| (\mathbf{A} \mathbf{x} - \mathbf{b})_{+} \|^2$，并不能直接将其分解成$p$个独立的子问题来分别求解，因为它们相互之间在目标函数上存在耦合。

基于原问题的特殊结构，作者提出了一种全新的块坐标下降方法。该方法通过将决策变量$\mathbf{x}$分为多个子块，通过按顺序循环迭代$\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_p$来最小化函数，从而逼近问题的最优解:
$$\mathbf{x}^{t+1}_j := \min_{\mathbf{x}_j \in \mathcal{X}_j} L(\mathbf{x}_j | \mathbf{x}^{t+1}_1, \ldots, \mathbf{x}^{t+1}_{j-1}, \mathbf{x}^t_{j+1}, \ldots, \mathbf{x}^t_p, \mathbf{\lambda}^k, {\rho}^k)$$
其中，其他块全部固定，只优化当前块$\mathbf{x}_j$。

论文中讨论了两种具体的块坐标更新方法：

- 近线性更新（Proximal Linear Update）
- 经典更新（Classical Update）

#### 4.1 近线性更新

给定固定的参数$\mathbf{\lambda} \geq \mathbf{0}$和$\rho > 0$，我们也可以计算$L(\mathbf{x}, \mathbf{\lambda}^k, {\rho}^k)$在$\mathbf{x}^t$处的梯度：
$$
g_j(\mathbf{x}^t) := \nabla_{\mathbf{x}_j}  L(\mathbf{x}, \mathbf{\lambda}, {\rho}) = \mathbf{c}_j + \mathbf{A}_j^{\intercal} \mathbf{\lambda} + \rho \mathbf{A}_j^{\intercal} (\mathbf{A} \mathbf{x}^t - \mathbf{b})_{+}
$$

近线性（Proximal Linear）方法可以利用这个一阶梯度$g_j(\mathbf{x}^t)$，将增广拉格朗日函数在点$\mathbf{x}^t$处线性近似为${g_j(\mathbf{x}^t)}^{\intercal} (\mathbf{x}_j - \mathbf{x}^t_j)$，同时忽略其高阶项。这种近似减少了计算的复杂性，同时保留了函数的主要变化趋势。接着，在线性化的目标函数中加入一个二次正则化项$\| \mathbf{x}_j - \mathbf{x}^t_j \|^2$，以促进解的稳定性，并防止更新步骤过大。可得：
$$
\mathbf{x}^{t+1}_j \in \arg \min_{\mathbf{x}_j \in \mathcal{X}_j} {g_j(\mathbf{x}^t)}^{\intercal} (\mathbf{x}_j - \mathbf{x}^t_j) + \frac{1}{\tau} \| \mathbf{x}_j - \mathbf{x}^t_j \|^2
$$
其中，$\tau$是步长参数，它是一个正数，用于控制更新的幅度；$\mathbf{x}^t_j$是在第$t$次迭代中块$j$的当前值。

该问题也可被写作：
$$
\mathbf{x}^{t+1}_j \in \arg \min_{\mathbf{x}_j \in \mathcal{X}_j} \| \mathbf{x}_j - (\mathbf{x}^t_j - \tau g_j(\mathbf{x}^t) ) \|^2
$$

一般来说，由于目标函数中的二次项，这个子问题并不容易求解。得益于$\mathbf{x}$是二元（0-1）变量，我们有
$$\| \mathbf{x}_j \|^2 = \mathbf{1}^{\intercal} \mathbf{x}_j,
$$
这里的$\mathbf{1}$是一个所有元素都是1的向量。这允许我们将原问题的二次部分线性化，进一步简化为一个整数线性问题，从而简化了求解过程：
$$
\begin{align}
\mathbf{x}^{t+1}_j & \in \arg \min_{\mathbf{x}_j \in \mathcal{X}_j} \mathbf{1}^{\intercal} \mathbf{x}_j - 2 {\mathbf{x}_j^t}^{\intercal} \mathbf{x}_j + 2 \tau {g_j(\mathbf{x}^t)}^{\intercal} \mathbf{x}_j \\
& = \arg \min_{\mathbf{x}_j \in \mathcal{X}_j} { \Big( \tau {g_j(\mathbf{x}^t)}  + \frac{\mathbf{1}}{2} - \mathbf{x}_j^t \Big)}^{\intercal} \mathbf{x}_j
\end{align}
$$

此外，步长$\tau$的选择对于算法的收敛性和效率至关重要。如果步长过大或过小，都可能影响算法的性能。

#### 4.2 经典更新

在论文中，为了简化子问题的求解，作者提出了一组特定的假设条件：系数矩阵 $\mathbf{A}$ 的元素仅包含 $0$ 或 $1$，同时右侧向量 $\mathbf{b} = \mathbf{1}$。虽然这些假设条件相当严格，但它们在许多优化问题中十分常见，如一些指派问题、路径问题和调度问题等。

在这些假设下，对于块坐标下降过程中的每一块 $\mathbf{x}_j$，不等式 $\mathbf{A}_j \mathbf{x}^t_j \leq \mathbf{1}$ 始终有效。利用这一特殊结构，作者能够推导出子问题的简化形式，将原本复杂的子问题线性化，并确保了算法在每一步都能产生满足约束条件的解：
$$
\mathbf{x}^{t+1}_j \in \arg \min_{\mathbf{x}_j \in \mathcal{X}_j} {\Big( \mathbf{c}_j + \mathbf{A}_j^{\intercal} \mathbf{\lambda} + \rho \mathbf{A}_j^{\intercal} (\sum_{l \neq j}^p \mathbf{A}_l \mathbf{x}_l^t - \frac{\mathbf{1}}{2})_+ \Big)}^{\intercal} \mathbf{x}_j
$$

通过这种简化，增广拉格朗日函数 $\min_{\mathbf{x} \in \mathcal{X}} L(\mathbf{x}, \mathbf{\lambda}, {\rho})$ 被分解为一系列具有线性目标函数的子问题。这种分解不仅降低了问题的复杂性，而且使得每个子问题更容易求解，从而提高了整个算法的效率。

#### 4.3 步长和收敛性

收敛性分析部分提供了理论保证，表明在适当的步长和初始条件下，块坐标下降法能够收敛到问题的最优解或者满足一定精度要求的可行解。

首先，如果起始点满足$\mathbf{x}^0 \in \mathcal{X}$，则该方法始终是可执行的，并在有限次迭代后终止，得到增广拉格朗日松弛问题$d(\mathbf{\lambda}, \rho)$的块最优解，即当保持其他块的变量不变时，该块的变量选择是最优的。

文中证明了证明了增广拉格朗日函数的梯度是Lipschitz连续的，这是全局收敛性的关键。对于全部$\mathbf{x}, \mathbf{x}^{\prime} \in \mathcal{X}$，$\kappa = \rho {\| \mathbf{A} \|}^2$有:
$$
\| \nabla L(\mathbf{x}, \mathbf{\lambda}, {\rho}) - L(\mathbf{x}^{\prime}, \mathbf{\lambda}, {\rho}) \| \leq \kappa \| \mathbf{x} - \mathbf{x}^{\prime} \|
$$

论文中进一步指出，近线性更新中步长 $\tau$ 需要要根据Lipschitz常数$\kappa$来确定。适当的步长选择确保了每次迭代后拉格朗日增广函数值的非增性。如果步长过大，算法可能无法收敛；步长过小，则迭代点 $\mathbf{x}^t$ 在近线性更新后可能没有变化。论文提出了一个理论界限 $0 < \tau < \frac{1}{2 \kappa}$，以确保算法能够在有限的迭代次数内收敛到 $\tau$-静止点，并证明了在足够大的步长下，$\tau$-静止点也是局部最小值，从而有助于确保块坐标下降算法的全局收敛性。

### 5 在块坐标下降法中寻找可行解

块坐标下降是一种通过逐步优化每个子块来解决大规模优化问题的迭代方法。但由于全局约束$\mathbf{A} \mathbf{x} \leq \mathbf{b}$的存在，子块之间存在相互依赖，单独优化某个子块时会导致全局解的不可行。尽管每次迭代可能会产生一个不可行的解，但算法可以结合额外的策略来修正这些不可行的解。这些策略不仅有助于在迭代过程中逐步确保解的可行性，还可以提供问题的上界估计，从而帮助计算对偶间隙。

因此，作者还讨论了在块坐标下降法寻找高质量全局可行解的两种技术：

- 扫掠技术（Sweeping Technique）
- 打包技术（Packing Technique）

这两种技术都依赖于为每个块$j$构造一个“解池”$V_j^k = \{ \mathbf{x}_j^1, \mathbf{x}_j^2, \ldots, \mathbf{x}_j^k \}$。解池收集了在先前迭代中获得的所有可行解，提供了丰富的备选解集。

#### 5.1 扫掠技术

扫掠技术是一种直观且高效的方法，用于在块坐标下降法中寻找可行解。它通过顺序检查每个块的解，并进行选择或重置，以确保全局解的可行性。

对于优化问题中的每个块 $j$，该技术按照以下步骤操作：

1. 从解池 $V_j^k$ 中选择一个候选解$\mathbf{v}_j$。
2. 检查当前候选解$\mathbf{v}_j$是否满足问题的全局约束$\mathbf{A}_j \mathbf{v}_j + \sum_{l \neq j}^p \mathbf{A}_l \hat{\mathbf{x}}^k_l \leq \mathbf{b}$。
3. 如果候选解$\mathbf{v}_j$是可行的，则将其保留作为该块的当前解$\hat{\mathbf{x}}^k_j := \mathbf{v}_j$。
4. 如果候选解$\mathbf{v}_j$不可行，将该块的解临时设为零（或另一个可行的初始状态），并继续检查序列中的下一个块。

在大规模问题或复杂约束条件下，扫掠技术可能无法找到质量较高的解。

#### 5.2 打包技术

打包技术用于从每个块的候选解集合 ${V_1^k, V_2^k, \ldots, V_p^k}$ 中选择出一组解 ${\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_p}$，这组解满足全局约束，共同构成问题的最优解。

对于每个模块 $j$，我们定义 $X_j^k := [ \mathbf{x}^1_j; \mathbf{x}^2_j; \ldots; \mathbf{x}^k_j]$ 为一个矩阵，其中包含了前 $k$ 次迭代中生成的所有解。接着，引入一组二值选择变量 $\mathbf{\mu}_j \in {\{ 0, 1 \}}^k$，使得 $\hat{\mathbf{x}}_j = X_j^k \mathbf{\mu}_j$。这里，二进制变量 $\mathbf{\mu}_j$ 决定了选择哪个解。

这种方法将问题转化为一个受限的主问题，其中块坐标下降法过程中生成的解通过列生成技术逐步加入到主问题中，以寻找最佳的解组合。优化问题的数学表达如下：
$$
\begin{align}
\min & \sum_{j=1}^{p} \mathbf{c}_j^{\intercal} X_j^k \mathbf{\mu}_j \\
\text{s.t.} \quad & \mathbf{\mu}_j^{\intercal} \mathbf{1} \leq 1, \quad \forall j = 1, \dots, p \\
& \sum_{j=1}^{p} \mathbf{A}_j X_j^k \mathbf{\mu}_j \leq \mathbf{b}
\end{align}
$$

由于全局约束（10）的存在，候选解之间存在冲突（即它们的组合无法满足全局约束），这些冲突使得这些解无法同时作为最终解的一部分。因此，可以将问题进一步转化为一个最大独立集（Maximal Independent Set, MIS）问题。最大独立集是指在一个图中找到一个不相邻节点的最大子集，这些节点之间没有边相连。

该方法需要构建一个冲突图$F = (V, E)$，其中节点集$V$表示所有候选解，而边集$E$表示这些候选解之间的冲突关系。如果两个解之间存在冲突（即它们无法同时满足某些约束），则在它们之间建立一条边。

然而，最大独立集问题仍然是一个NP难问题。文中提到实际上求解的是松弛后的最大独立集问题。通过松弛，可以在一定程度上简化计算，并找到一个可行解。但是文中**没有**明确指出采用了哪种具体的算法来处理这一松弛问题。


### 实验结果和总结

在论文的实验部分，作者对所提出的增广拉格朗日方法进行了数值实验，以验证其在求解块结构整数规划问题中的有效性。实验使用了几个具有实际应用背景的场景，具体包括火车时刻表安排问题和车辆路径问题。这些问题其规模大且包含整数变量。论文展示了增广拉格朗日算法与其他现有方法相比在求解效率和解的质量上的优势。

- 求解效率：该算法在多个实例中显示出快速的收敛速度，能够在较少的迭代次数内达到满意的解。
- 解的质量：该算法能够找到高质量的可行解，并且在许多情况下接近全局最优解。
- 稳定性：在不同的测试案例中，该算法表现出良好的稳定性和鲁棒性，即使在问题规模较大或数据波动时也能保持性能。


#### 数据集

在论文的 "5.1 Capacitated Vehicle Routing Problem" 部分，作者提到了使用 Solomon 数据集来评估他们提出的算法的性能。Solomon 数据集是一组广泛用于评估车队路径问题（Vehicle Routing Problem, VRP）算法的标准测试实例。这些数据集包含了不同的参数，例如车辆数量、客户数量、服务时间窗口以及车辆容量限制等。这里，我们使用Solomon 数据集的C102实例进行算法的复现。

这次复现工作选取了Solomon数据集中的C102实例。C102实例包含100个客户（不包括始发点）。在论文中，作者考虑了C102实例的不同变体，以探究算法在不同规模问题上的表现：

- C102-25：选择了前25个客户，分配了3辆车辆。
- C102-50：选择了前50个客户，同样分配了3辆车辆。
- C102-100：选择了全部100个客户，分配了10辆车辆。

数据集存储于文件`c102.txt`中，并由`data.py`脚本进行读取和处理。

**关于车辆的最大载货量，原文中没有明确提及具体设置。根据作者的实验结果，我们可以暂时将车辆的最大载货量设置为$200$单位。**

#### 数学模型

对于带容量限制的车辆路径问题（Capacitated Vehicle Routing Problem, CVRP），我们定义以下符号：

- $V$: 节点集合，其中$0$代表仓库，$1$到$n$代表客户。
- $E$：边集合，连接节点。
- $\mathbb{N}_p$: 车辆的编号集合
- $x^j_{st}$: 若车辆$j$行驶路径 $(s, t)$被选中则为1，否则为0。
- $w^j_s$: 车辆$j$开始在客户$s$服务的时间。
- $d_{st}$: 节点$s$到节点$t$的距离。
- $c_s$: 客户$s$的需求量。
- $C$: 车辆的最大载货量。
- $T_{st}$: 从节点$s$到节点$t$需要的时间。
- $a_s$ 和 $b_s$: 客户$s$的服务时间窗口的开始和结束时间。
- $M$: 一个足够大的常数，用于确保时间窗口约束。**文中没有提及$M$的设置，我们这里将$M$选择为服务时间窗口的结束时间的最大值。**

有数学模型：
$$
\begin{align}
\min & \sum_{j \in \mathbb{N}_p} \sum_{(s,t) \in E} d_{st} x^j_{st} \\
\text{s.t.} \quad & \sum_{j \in \mathbb{N}_p} \sum_{t \in V : t \neq s} x^j_{st} = 1, \quad & \forall s \in V \setminus \{0\} \\
& \sum_{t \in V \setminus \{s\}} x^j_{st} = \sum_{t \in V \setminus \{s\}} x^j_{ts}, \quad & \forall s \in V, \forall j \in \mathbb{N}_p \\
& \sum_{t \in V \setminus \{0\}} x^j_{0t} = 1, \quad & \forall j \in \mathbb{N}_p \\
& \sum_{s \in V} \sum_{t \in V \setminus \{s\}} c_s x^j_{st} \leq C, \quad & \forall j \in \mathbb{N}_p \\
& w^j_s + T_{st} - M(1 - x^j_{st}) \leq w^j_t, \quad & \forall (s, t) \in E, \forall j \in \mathbb{N}_p \\
& a_s \leq w^j_s \leq b_s, \quad & \forall s \in V, \forall j \in \mathbb{N}_p \\
& x^j_{st} \in \{0, 1\} \quad & \forall (s, t) \in E, \forall j \in \mathbb{N}_p, \\
& w^j_s \geq 0 \quad & \forall s \in V, \forall j \in \mathbb{N}_p \\
\end{align}
$$
其中，目标函数（11）最小化所有车辆路径的总距离；约束条件（12）要求每个客户恰被一辆车服务一次；约束条件（13）保证每个车辆在每个节点的进入和离开流量相等；约束条件（14）确保每个车辆从仓库出发并最终返回；约束条件（15）令每辆车的载货量不超过其容量；约束条件（16）限制车辆在客户节点的服务开始时间必须在时间窗口内；约束条件（17）规定车辆到达每个客户节点的时间必须在该节点的时间窗口内。

**注意，尽管原文中没有提及，此处约束条件（16）中$s$不能为始发站，否则会导致模型无解。**

此Gurobi模型见`gurobi.py`，通过求解器Gurobi求解，作为算法性能的基准测试。求解结果与论文中的结果一致，验证了模型的准确性。

#### 增广拉格朗日函数

在这里，我们发现只有约束条件（12）是全局约束。通过松弛约束（12），问题将被分解为每辆车单独的路径子问题。因此，我们可以实施增广拉格朗日方法。

假设 $\mathbf{\lambda}$ 是约束（12）的拉格朗日乘子，$\rho > 0$ 是增广项的系数，则增广拉格朗日函数可以表示为：
$$
f(\mathbf{x}) + \sum_{s \in V \setminus \{0\}} {\lambda}_s \left(\sum_{j \in \mathbb{N}_p} \sum_{t \in V : t \neq s} x^j_{st} - 1 \right) + \frac{\rho}{2} \sum_{s \in V \setminus \{0\}} \left(\sum_{j \in \mathbb{N}_p} \sum_{t \in V : t \neq s} x^j_{st} - 1 \right)^2
$$
其中$f(\mathbf{x})$是原目标函数。

全局约束的右侧值 $\mathbf{b} = 1$。对于每辆车 $j$，约束矩阵块 $\mathbf{A}_j$ 是一个 $(|V|-1) \times |E|$ 的矩阵。具体来说，$\mathbf{A}j$ 的行对应于每个客户节点 $s$（不包括仓库），而列对应于每条可能的路径 $x^j_{st}$。如果路径 $(s, t)$ 被包含在约束中，则对应的矩阵元素为 $1$，否则为 $0$。

假设有三个客户（节点1、节点2、节点3）以及一个仓库（节点0），矩阵$\mathbf{A}_j$可以表示为：
$$
\mathbf{A}_j =
\begin{pmatrix}
0 & 0 & 0 & 1 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 1
\end{pmatrix}
$$

实现增广拉格朗日方法的代码可以在`alm.py`中找到，而块坐标下降法的实现则在`bcd.py`中。

#### 初始值

在文章中，作者并未详细讨论车辆路径的增广拉格朗日方法中的初始解$\mathbf{x}^0$、拉格朗日乘子初始值$\mathbf{\lambda}^0$和惩罚参数初始值$\mathbf{\lambda}^0$。

为了寻找一个合适的初始解 $\mathbf{x}^0$，我们采用了一种贪心算法（实现代码见`greedy.py`）。该算法的核心思想是：对于当前节点 $s$，在确保时间窗约束和车辆容量限制的前提下，选择下一个节点 $t$，使得 $b_t - d_{st}$ 的值最大。这里，$b_t$ 表示节点 $t$ 的时间窗约束的结束时间，而 $d_{st}$ 是从节点 $s$ 到节点 $t$ 的距离（时间）。通过这种方式，我们能够为每辆车构建一个初始的路径。通过这种方法，我们能够为每辆车构建一个初始的路径。虽然这个算法不能保证满足所有全局约束，但它确实满足了块约束 $\mathbf{x}_j \in \mathcal{X}_j$，符合文章中收敛性的要求。

初始拉格朗日乘子 $\mathbf{\lambda}^0$设置为零向量；初始惩罚参数$\mathbf{\lambda}^0$设置为1。

#### 结论

增广拉格朗日方法为解决大规模车辆提供了一种有效的数学规划框架，通过合理选择初始值和调整算法参数，可以显著提高求解效率和解的质量。在实际应用中，增广拉格朗日方法展现出了与论文中提出的结果相一致的性能：它能够比其他现有方法更快地收敛，尤其在处理大规模问题时，这一优势更为明显。


