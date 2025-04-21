# 二维FSI程序求解



## Sec 1. 交错网格上分布的物理信息

<img src="C:\Users\wangz\AppData\Roaming\Typora\typora-user-images\image-20250414101936980.png" alt="image-20250414101936980" style="zoom:50%;" />



## Sec 2. N-S方程的求解方法：一阶投影法

将时间推进分解成三个子步长，压力不通过时间推进求解：

Step 1：预算步：
$$
\begin{array}{l}
\frac{\mathbf{V}^{*}-\mathbf{V}^{n}}{\Delta t}+\left(\mathbf{V} \cdot \nabla \mathbf{V}-\frac{1}{\operatorname{Re}} \nabla^{2} \mathbf{V}\right)^{n}=0
\end{array}
$$
Step 2：压力修正步
$$
\begin{array}{l}
\frac{\mathbf{V}^{n+1}-\mathbf{V}^{*}}{\Delta t}+\nabla p=0 \to   \nabla^{2} p=\frac{1}{\Delta t} \nabla \cdot \mathbf{V}^{*} \quad   \\
\nabla \cdot \mathbf{V}^{n+1}=0
\end{array}
$$
Step 3：最终步

$$
\frac{\mathbf{V}^{n+1}-\mathbf{V}^{*}}{\Delta t}+\nabla p=0 \quad \text { 得到 } \mathbf{n}+\mathbf{1} \text { 时刻的 } \mathbf{V}
$$


## Sec 3. 交错网格上的动量方程离散（预算步）

(i，j)以压力网格为模板，u方向的动量离散公式：
$$
(u_{i+\frac{1}{2},j}^* - u_{i+\frac{1}{2},j}^n)/ \Delta t =  - u_{i+\frac{1}{2},j}^n \frac{\partial u }{\partial x}   - v_{i+\frac{1}{2},j}^n \frac{\partial u }{\partial y} + \frac{1}{\operatorname{Re}} \nabla^{2} u_{i+\frac{1}{2},j}^n
$$
其中，正方向的对流项的离散方案如下：
$$
u_{i+\frac{1}{2},j}^{} \frac{\partial u }{\partial x} = u_{i+\frac{1}{2},j}^{+} \frac{\partial u }{\partial x} + u_{i+\frac{1}{2},j}^{-} \frac{\partial u }{\partial x}，u_{i+\frac{1}{2},j}^{+，-} = \frac{u+,- |u|}{2},\ \frac{\partial u }{\partial x} = \frac{u_{i+\frac{3}{2},j} - u_{i+\frac{1}{2},j}}{\Delta x}
$$
而，对于剪切方向的离散方案，需要将v速度插值到u的网格模板上：
$$
v_{i+\frac{1}{2},j}^{} = 0.5 \times (v_{i+1,j}^{}+v_{i,j}^{})  = 0.5 \times ( \ 0.5 \times (v_{i+1,j-\frac{1}{2}}^{} + v_{i+1,j+\frac{1}{2}}^{})  \ + \ 0.5 \times (v_{i,j-\frac{1}{2}}^{} + v_{i,j+\frac{1}{2}}^{}))
$$


(i，j)以压力网格为模板，v方向的动量离散公式：
$$
(v_{i,j+\frac{1}{2}}^* - v_{i,j+\frac{1}{2}}^n)/ \Delta t =  - u_{i,j+\frac{1}{2}}^n \frac{\partial v }{\partial x}   - v_{i,j+\frac{1}{2}}^n \frac{\partial v }{\partial y} + \frac{1}{\operatorname{Re}} \nabla^{2} v_{i,j+\frac{1}{2}}^n
$$
其中，正方向的对流项的离散方案如下：
$$
v_{i,j+\frac{1}{2}}^{} \frac{\partial v }{\partial y} = v_{i,j+\frac{1}{2}}^{+} \frac{\partial v }{\partial y} + v_{i,j+\frac{1}{2}}^{-} \frac{\partial v }{\partial y}，v_{i,j+\frac{1}{2}}^{+，-} = \frac{v+,- |v|}{2},\ \frac{\partial v }{\partial y} = \frac{v_{i,j+\frac{3}{2}} - v_{i,j+\frac{1}{2}}}{\Delta y}
$$
同样，对于剪切方向的离散方案，需要将u速度插值到v的网格模板上：
$$
u_{i,j+\frac{1}{2}}^{} = 0.5 \times (u_{i,j+1}^{}+u_{i,j}^{})  = 0.5 \times ( \ 0.5 \times (u_{i+\frac{1}{2},j+1}^{} + u_{i-\frac{1}{2},j+1}^{})  \ + 0.5 \times (u_{i+\frac{1}{2},j}^{} + u_{i-\frac{1}{2},j}^{}))
$$


## Sec 4. 交错网格上的压力泊松方程离散（中间步）

压力泊松方程：
$$
\begin{array}{l}
\frac{\mathbf{V}^{n+1}-\mathbf{V}^{*}}{\Delta t}+\nabla p=0 \\
\nabla \cdot \mathbf{V}^{n+1}=0
\end{array}
$$
以压力网格为基准，交错网格上的离散结果：
$$
\begin{array}{l}
u_{i+1 / 2, j}^{n+1}-u_{i+1 / 2, j}^{*}+\Delta t / \Delta x\left(p_{i+1, j}-p_{i, j}\right)=0 \\
v_{i, j+1 / 2}^{n+1}-v_{i, j+1 / 2}^{*}+\Delta t / \Delta y\left(p_{i, j+1}^{*}-p_{i, j+1}\right)=0 \\
\left(u_{i+1 / 2, j}^{n+1}-u_{i-1 / 2, j}^{n+1}\right) / \Delta x+\left(v_{i, j+1 / 2}^{n+1}-v_{i, j-1 / 2}^{n+1}\right) / \Delta y=0
\end{array}
$$
将前两个式子带入到第三个方程（即不可压条件中），可以得到交错网格上，压力满足的泊松方程的离散形式：
$$
\frac{1}{\Delta t} ( \ \frac{u_{i+\frac{1}{2},j}^* - u_{i-\frac{1}{2},j}^*}{\Delta x}  \  +  \ \frac{v_{i,j+\frac{1}{2}}^* - v_{i,j-\frac{1}{2}}^*}{\Delta x} \ ) = \frac{p_{i+1,j} + p_{i-1,j} - 2p_{i,j} }{(\Delta x)^2}  + \frac{p_{i,j+1} + p_{i,j-1} - 2p_{i,j} }{(\Delta y)^2}
$$


## Sec 5. 交错网格上的压力矫正离散（最后步）

通过压力矫正可以获得速度场：
$$
\frac{\mathbf{V}^{n+1}-\mathbf{V}^{*}}{\Delta t}+\nabla p=0
$$
那么，交错网格上速度场u的矫正为：
$$
\frac{u_{i+\frac{1}{2},j}^{n+1}-u_{i+\frac{1}{2},j}^{*}}{\Delta t}+\frac{p_{i+1,j}-p_{i,j}}{\Delta x}=0
$$
同理，交错网格上速度场v的矫正为：
$$
\frac{v_{i,j+\frac{1}{2}}^{n+1}-v_{i,j+\frac{1}{2}}^{*}}{\Delta t}+\frac{p_{i,j+1}-p_{i,j}}{\Delta y}=0
$$

## Sec 6. 一些其他需要注意的地方！

1. 边界上的差分如何处理？

   示意图：| (n1) - - (n2) - - (n3)，这是靠近边界内部的三个节点，其中，n1，2，3代表节点，这里的n1就在边界上，先假设其网格节点是均匀的，都是dx，假设靠近壁面的速度可以用如下公式拟合（三个点用三次样条曲线）：
   $$
   u = a+by+cy^2 \\
   u_1 = a \\
   u_2 = a + b(dx) \\
   u_3 = a + b(2dx) + c(2dx)^2 \\
   b = \frac{-3u_1+4u_2-u_3}{2dx} \\
   c = \frac{u_3 + u_1 - 2u_2}{2dx^2} \\
   \frac{\partial u}{\partial y} = b+2cy, \frac{\partial u}{\partial y}|_{y=0} =b =  \frac{-3u_1+4u_2-u_3}{2dx}\\
   \frac{\partial^2 u}{\partial y^2} = 2c, \frac{\partial^2 u}{\partial y^2}|_{y=0} =2c =  \frac{u_3 + u_1 - 2u_2}{dx^2} \\
   $$
