
#include "LinearSystem.h"  

int main() {
    int nVar = 2;
    long nCell = 3;
    long cellStart = 0;

    std::vector<std::vector<long>> neighbors = {
        {0, 1},
        {0, 1, 2},
        {1, 2}
    };

    BlockCSRMatrix mat(nCell, cellStart, nVar, neighbors);

    // 创建一个 2x2 的块矩阵
    double block_data[] = {
        1.0, 2.0,
        3.0, 4.0
    };

    // 设置块 (0,1)
    mat.setBlock(0, 1, block_data);

    // 添加块 (1,2)，系数为 0.5
    mat.addBlock(1, 2, block_data, 0.5);

    // 打印单个块
    mat.printBlock(0, 1);
    mat.printBlock(1, 2);

    // 打印所有块
    mat.printAllBlocks();

    // 打印稀疏结构
    mat.printPattern();

    // 写入文件
    mat.printAllBlocksToFile("block_matrix.txt");

    return 0;
}

