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

## A time-dependent heat equation problem 

我们求解一个时间相关的二维热传导方程，具体的方程如下：
$$
\frac{\partial u}{\partial t} = D_0 \nabla^2 u + f(x, y).
$$
其中$u(t, x, y)  \ for  \ t \in [0, T] , (x, y) \in S = (0, 1) \times (0, 1),D_0 > 0, u(0, x, y) = 0, f(x, y) = 3 e^{-25 (x - 0.6)^2} \sin(2 \pi y)$.关于边界条件：水平方向上（$x = 0,x=1$）采用了Neumann条件，分别是：
$$
-\frac{\partial u}{ \partial x}= \gamma(y) = \sin(6\pi y) , x = 0\\
\frac{\partial u}{ \partial x}= 0,x = 1
$$
而在顶部和底部所采用的是，周期边界条件(上下底面的热量流动时平衡的):
$$
u(t,x,1) = u(t,x,0)
$$
关于这样的热传导方程，我们有一个定性的分析是：
$$
\int \frac{\partial u}{\partial t}  \ d \Omega = \int D_0 \nabla^2 u + f(x, y)  \ d \Omega \\
 = \int D_0 \nabla \cdot (\nabla u) + f(x, y)  \ d \Omega \\
  = \int D_0  \nabla u  \cdot n \ d S + \int  f(x, y)  \ d \Omega  \\
   = \int D_0  \frac{\partial u}{ \partial x}   \ d y + \int  f(x, y)  \ d \Omega \\ = 
   \int -  D_0  \sin(6\pi y)   \ d y + \int  3 e^{-25 (x - 0.6)^2} \sin(2 \pi y)  \ d \Omega \\
    = 0 = \frac{\partial }{\partial t}\int u\ d \Omega
$$
说明热能守恒，所以我们希望在数值离散的时候，也能保证这一种特性！

使用二阶中心差分来离散laplacian项，则有离散后的方程为：
$$
u^{'}_{i,j} = D_0 \frac{u_{i-1,j}  \ -  \ 2u_{i,j} \ + \ u_{i+1,j}}{h_{x}^{2}}  + D_0 \frac{u_{i,j-1}  \ -  \ 2u_{i,j} \ + \ u_{i,j+1}}{h_{y}^{2}} + f_{i,j} \\
  which \ \ is  \ \ equal \ \ to \  :\ u' = g(u)
$$
