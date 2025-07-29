#!/bin/bash

echo "==> 编译程序..."
mpicc ij_matrix_example.c \
  -I/mnt/d/MyProject/HYPRE/install/include \
  -L/mnt/d/MyProject/HYPRE/install/lib \
  -lHYPRE -lm -o ij_matrix_example

if [ $? -ne 0 ]; then
  echo "❌ 编译失败"
  exit 1
fi

echo "==> 运行程序..."
mpirun -n 1 ./ij_matrix_example > ij_matrix_example.log 2>&1

echo "✅ 输出已保存到 ij_matrix_example.log"




