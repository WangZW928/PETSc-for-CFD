# FDM for Numerical Computation by PETSc

## 2D均匀网格 for Coupled reaction-diffusion equations 

耦合反应-扩散方程常用于描述多个物理量或物种（如化学物质、物理场或生物种群等）之间的相互作用和扩散过程。它们通常包括多个方程，每个方程描述了一个物种或物理量的变化，这些物种或物理量通过反应项和扩散项相互耦合。其一般的形式可以写为：

$$
\frac{\partial u}{\partial t} = L u + X(u), \quad u(0) = f,
$$
其中$u = u(t, x) \in \mathbb{R}^k$，$X $是定义在$\mathbb{R}^k$的实向量场，且$L$是一个二阶微分算子。在数值计算时我们考虑在欧式空间中的如下方程：
$$
\frac{\partial u}{\partial t} &= D_u \nabla^2 u - uv^2 + \phi(1 - u), \\
\frac{\partial v}{\partial t} &= D_v \nabla^2 v + uv^2 - (\phi + \kappa) v.
$$
其中$u(t, x, y)$，$v(t, x, y)$是化学反应物浓度，两者在空间中都发生了耗散，且耗散系数分别为：$D_u > 0，D_v > 0$。反应物$u$会生成$v$，生成速率为$uv^2$，另外反应物$u$还有额外的补充：$\phi > 0$，而反应物$v$会生成额外的组分。计算域为：$S = [0, L) \times [0, L) $。
>**注意了，这种时间演化的非线性方程一定要注意解是否能稳定（收敛）**
