# 编译器设置
CXX := mpic++
CXXFLAGS := -std=c++17 -O2 -Wall -Iinclude -MMD -MP

# 目录
SRC_DIR := src
OBJ_DIR := build
BIN := main

# 自动收集源文件 / 对象文件 / 依赖文件
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))
DEPS := $(OBJS:.o=.d)

# 提取目标名，如 main.cpp -> main
TARGETS := $(notdir $(basename $(SRCS)))

# 默认目标：编译所有并链接主程序
all: $(BIN)

# 链接最终可执行文件
$(BIN): $(OBJS)
	$(CXX) $^ -o $@

# 单独编译文件，例如 make util 生成 ./util 可执行文件
$(TARGETS): %: $(OBJ_DIR)/%.o
	$(CXX) $< -o $@

# 编译规则（源文件 -> 对象文件）
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 创建目标目录
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# 清理所有编译中间产物和可执行文件
clean:
	rm -rf $(OBJ_DIR) *.o *.d $(BIN) $(TARGETS)

# 伪目标
.PHONY: all clean $(TARGETS)

# 自动包含头文件依赖
-include $(DEPS)
