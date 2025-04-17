# 设置 PETSc 相关环境变量
export PETSC_DIR=/home/wangzw/petsc
export PETSC_ARCH=arch-linux-c-debug

# 获取 PETSc 库路径
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# 添加 C++ 编译器选项
CXXFLAGS += -pedantic -std=c++11  # 使用 C++11 标准（你也可以使用 c++14 或 c++17）

# 如果你想在 C++ 中使用 MPI，使用 mpicxx 或 g++
CXX := mpicxx  # 或者使用 g++

# 自动获取当前目录下所有的 .cpp 文件
SRCS := $(wildcard *.cpp)

# 将所有 .cpp 文件编译为 .o 文件
OBJS := $(SRCS:.cpp=.o)

# 最终的可执行文件名称
TARGET := test

# 默认目标（编译可执行文件）
all: $(TARGET)

# 生成可执行文件
$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) ${PETSC_LIB}


# 完全清理（包括所有临时文件）
distclean: 
	rm -f *~ 

# 标记目标，不是文件名
.PHONY: distclean  all
