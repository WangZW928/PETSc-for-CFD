
# 设置 PETSc 相关环境变量
export PETSC_DIR=/home/wangzewei/PETSc/petsc
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# 添加 C++ 编译器选项
CXXFLAGS += -pedantic -std=c++11  # 使用 C++11 标准（你也可以使用 c++14 或 c++17）

# 如果你想在 C++ 中使用 MPI，使用 mpicxx 或 g++
CXX := mpicxx  # 或者使用 g++

# 指定源文件扩展名为 .cpp
SRCS := $(wildcard *.cpp)

# 编译后的对象文件
OBJS := $(SRCS:.cpp=.o)

# 目标：生成可执行文件
e: $(OBJS)
	-${CXX} -o e $(OBJS) ${PETSC_LIB}
	${RM} $(OBJS)

# 测试目标（可以根据需要调整）
rune_1:
	-@../testit.sh e "" 1 1

test_e: rune_1

test: test_e

# 清理目标
clean::
	@rm -f *~ e $(OBJS) *tmp

# 完全清理
distclean: clean

# 标记目标，不是文件名
.PHONY: clean distclean rune_1 test_e test

