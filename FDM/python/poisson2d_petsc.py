import numpy as np
from petsc4py import PETSc

# 创建上下文
comm = PETSc.COMM_WORLD
rank = comm.getRank()

# 设置网格尺寸（n x n）
n = 32
h = 1.0 / (n + 1)

# 创建矩阵 A
A = PETSc.Mat().createAIJ([n * n, n * n], nnz=5, comm=comm)
A.setFromOptions()
A.setUp()

# 创建右端项 b 和解向量 x
b = PETSc.Vec().createSeq(n * n, comm=comm)
x = b.duplicate()

# 组装 A 和 b
for i in range(n):
    for j in range(n):
        row = i * n + j
        xi = (i + 1) * h
        yj = (j + 1) * h
        b[row] = h**2 * 2 * np.pi**2 * np.sin(np.pi * xi) * np.sin(np.pi * yj)

        A.setValue(row, row, 4.0)
        if i > 0:    A.setValue(row, row - n, -1.0)
        if i < n-1:  A.setValue(row, row + n, -1.0)
        if j > 0:    A.setValue(row, row - 1, -1.0)
        if j < n-1:  A.setValue(row, row + 1, -1.0)

A.assemble()
b.assemble()

# 解线性系统 Ax = b
ksp = PETSc.KSP().create(comm)
ksp.setOperators(A)
ksp.setType('cg')
ksp.getPC().setType('jacobi')
ksp.setFromOptions()
ksp.solve(b, x)

# 计算误差
error = 0.0
for i in range(n):
    for j in range(n):
        idx = i * n + j
        xi = (i + 1) * h
        yj = (j + 1) * h
        u_exact = np.sin(np.pi * xi) * np.sin(np.pi * yj)
        error += (x[idx] - u_exact) ** 2

error = np.sqrt(error / (n * n))
if rank == 0:
    print(f"✅ L2 误差 = {error:.3e}")
