#!/bin/bash

echo "==> 编译程序..."

mpicxx testMatrix.cpp \
  -I/mnt/d/MyProject/HYPRE/install/include \
  -L/mnt/d/MyProject/HYPRE/install/lib \
  -lHYPRE -lm -o test_matrix

if [ $? -ne 0 ]; then
  echo "❌ 编译失败"
  exit 1
fi

echo "==> 运行程序..."
mpirun -n 1 ./test_matrix > test_matrix.log 2>&1

echo "✅ 输出已保存到 test_matrix.log"
