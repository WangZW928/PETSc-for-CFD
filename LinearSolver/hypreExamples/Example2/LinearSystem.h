
#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

using mat_t = std::array<double, 25>;  // 支持 nVar <= 5 时的块
    
class BlockCSRMatrix {
public:
    int nVar = 0;
    long nCell = 0;
    long globalCellStart = 0;

    std::vector<long> row_ptr;
    std::vector<long> col_idx;
    std::vector<double> blocks;

    BlockCSRMatrix(long nCell_, long cellStart, int nVar_, const std::vector<std::vector<long>>& neighbors)
        : nVar(nVar_), nCell(nCell_), globalCellStart(cellStart) {
        row_ptr.resize(nCell + 1);
        long total_blocks = 0;

        for (long i = 0; i < nCell; i++) {
            row_ptr[i] = total_blocks;
            total_blocks += neighbors[i].size();
        }
        row_ptr[nCell] = total_blocks;

        col_idx.resize(total_blocks);
        long idx = 0;
        for (long i = 0; i < nCell; i++) {
            for (long j : neighbors[i]) {
                col_idx[idx++] = j;
            }
        }

        blocks.resize(total_blocks * nVar * nVar, 0.0);
    }

    // 查找结果封装
    struct BlockSearchResult {
        long start;
        long block_idx;
        long num_blocks;
        double* block_ptr;  // 指向块数据的指针
    };

    BlockSearchResult findBlock(long cell_i, long cell_j) {
        if (cell_i < globalCellStart || cell_i >= globalCellStart + nCell)
            throw std::out_of_range("cell_i out of range");

        long local_i = cell_i - globalCellStart;
        long start = row_ptr[local_i];
        long end = row_ptr[local_i + 1];
        long num_blocks = end - start;

        for (long i = start; i < end; i++) {
            if (col_idx[i] == cell_j) {
                long block_idx = i - start;
                double* base = blocks.data() + block_idx * nVar * nVar;
                return { start, block_idx, num_blocks, base };
            }
        }

        throw std::runtime_error("Block(" + std::to_string(cell_i) + ", " + std::to_string(cell_j) + ") not found.");
    }

    void setBlock(long cell_i, long cell_j, const double* data) {
        auto result = findBlock(cell_i, cell_j);

        for (int i = 0; i < nVar; ++i) {
            double* row_ptr = result.block_ptr + i * result.num_blocks * nVar + result.block_idx * nVar;
            for (int j = 0; j < nVar; ++j) {
                row_ptr[j] = data[i * nVar + j];
            }
        }
    }

    void addBlock(long cell_i, long cell_j, const double* data, double factor = 1.0) {
        auto result = findBlock(cell_i, cell_j);

        for (int i = 0; i < nVar; ++i) {
            double* row_ptr = result.block_ptr + i * result.num_blocks * nVar + result.block_idx * nVar;
            for (int j = 0; j < nVar; ++j) {
                row_ptr[j] += factor * data[i * nVar + j];
            }
        }
    }

    void printBlock(long cell_i, long cell_j) const {
        auto result = const_cast<BlockCSRMatrix*>(this)->findBlock(cell_i, cell_j);
        std::cout << "Block(" << cell_i << ", " << cell_j << "):\n";

        for (int i = 0; i < nVar; ++i) {
            const double* row_ptr = result.block_ptr + i * result.num_blocks * nVar + result.block_idx * nVar;
            for (int j = 0; j < nVar; ++j) {
                std::cout << row_ptr[j] << " ";
            }
            std::cout << "\n";
        }
    }

    void printAllBlocks() const {
        std::cout << "===== Matrix Blocks =====\n";
        for (long i = 0; i < nCell; ++i) {
            long start = row_ptr[i];
            long end = row_ptr[i + 1];
            long num_blocks = end - start;

            for (long k = 0; k < num_blocks; ++k) {
                long j = col_idx[start + k];
                std::cout << "Block(" << globalCellStart + i << ", " << j << "):\n";
                const double* base = blocks.data() + start * nVar * nVar;

                for (int row = 0; row < nVar; ++row) {
                    const double* row_ptr = base + row * num_blocks * nVar + k * nVar;
                    for (int col = 0; col < nVar; ++col) {
                        std::cout << row_ptr[col] << " ";
                    }
                    std::cout << "\n";
                }
                std::cout << "\n";
            }
        }
    }

    void printAllBlocksToFile(const std::string& filename) const {
        std::ofstream out(filename);
        if (!out) {
            std::cerr << "❌ Cannot open file: " << filename << "\n";
            return;
        }

        out << "===== Matrix Blocks =====\n";
        for (long i = 0; i < nCell; ++i) {
            long start = row_ptr[i];
            long end = row_ptr[i + 1];
            long num_blocks = end - start;

            for (long k = 0; k < num_blocks; ++k) {
                long j = col_idx[start + k];
                out << "Block(" << globalCellStart + i << ", " << j << "):\n";
                const double* base = blocks.data() + start * nVar * nVar;

                for (int row = 0; row < nVar; ++row) {
                    const double* row_ptr = base + row * num_blocks * nVar + k * nVar;
                    for (int col = 0; col < nVar; ++col) {
                        out << row_ptr[col] << " ";
                    }
                    out << "\n";
                }
                out << "\n";
            }
        }

        std::cout << "✅ Matrix written to: " << filename << "\n";
    }

    void printPattern(double threshold = 1e-12) const {
        long totalRows = nCell * nVar;
        long totalCols = totalRows;  // 简化：假设是对称结构

        std::vector<std::string> pattern(totalRows, std::string(totalCols, ' '));

        for (long i = 0; i < nCell; ++i) {
            long start = row_ptr[i];
            long end = row_ptr[i + 1];
            long num_blocks = end - start;

            for (long k = 0; k < num_blocks; ++k) {
                long j = col_idx[start + k];
                long row_offset = i * nVar;
                long col_offset = (j - globalCellStart) * nVar;

                const double* base = blocks.data() + start * nVar * nVar;

                for (int ii = 0; ii < nVar; ++ii) {
                    const double* row_ptr = base + ii * num_blocks * nVar + k * nVar;
                    for (int jj = 0; jj < nVar; ++jj) {
                        if (std::fabs(row_ptr[jj]) > threshold) {
                            pattern[row_offset + ii][col_offset + jj] = '*';
                        }
                    }
                }
            }
        }

        std::cout << "Matrix pattern (" << totalRows << " x " << totalCols << "):\n";
        for (const auto& line : pattern) {
            std::cout << line << "\n";
        }
    }

    void validate() const {
        if (row_ptr.size() != static_cast<size_t>(nCell + 1))
            throw std::runtime_error("Invalid row_ptr size");

        if (col_idx.size() != static_cast<size_t>(row_ptr[nCell]))
            throw std::runtime_error("col_idx size does not match row_ptr");

        if (blocks.size() != static_cast<size_t>(col_idx.size() * nVar * nVar))
            throw std::runtime_error("blocks size inconsistent with structure");
    }
};





/*

class LinearSystem
{
public:
    HYPRE IJMatrix A;
    HYPRE IJVector b;
    HYPRE IJVector x;   
    HYPRE_ParCSRMatrix parcsrA;
    HYPRE_ParCSRVector parB;
    HYPRE_ParCSRVector parX;
	LinearSystem();
	~LinearSystem();

    void addJacobiBlock(HYPRE_BigInt i, HYPRE_BigInt j,double** jacobian);
    void subJacobiBlock(HYPRE_BigInt i, HYPRE_BigInt j,double** jacobian);

private:
    HYPRE_BigInt iLower = 0, iUpper = 0;
    HYPRE_Int nRows = 0;
    std::vector<HYPRE_Int> nCols;
    std::vector<HYPRE_BigInt> rows;
    std::vector<HYPRE_BigInt> cols;
    HYPRE_BigInt globalCellStart = 0;

    HYPRE_Int nVar = 0;
    unsigned short nDim = 0;
    bool implicit = false;
    int myid = -1;

    std::unique_ptr<BlockCSRMatrix> matrix;
    std::vector<double> buffer;
    std::vector<double> bPtr,xPtr;

    // Rule 
    LinearSystem(const LinearSystem&) = delete;
    LinearSystem& operator=(const LinearSystem&) = delete;
    LinearSystem(LinearSystem&&) noexcept = default;
    LinearSystem& operator=(LinearSystem&&) noexcept = default;



};

inline LinearSystem::LinearSystem() {
    // Constructor implementation
}

inline LinearSystem::~LinearSystem() {
    // Destructor implementation
}

void LinearSystem::addJacobiBlock(HYPRE_BigInt i, HYPRE_BigInt j, mat_t jacobian) {
    // Implementation for adding a Jacobi block
    for (int ielem = 0; ielem < 25; ++ielem) {
        // Example of adding the block to the matrix
        buffer[ielem] = jacobian[ielem];
    }
    matrix->addBlock(i, j, buffer.data());

}

void LinearSystem::subJacobiBlock(HYPRE_BigInt i, HYPRE_BigInt j, mat_t jacobian) {
    // Implementation for subtracting a Jacobi block
    for (int ielem = 0; ielem < 25; ++ielem) {
        buffer[ielem] = jacobian[ielem];
    }
    matrix->addBlock(i, j, buffer.data(), -1.0);
}

*/


#endif // LINEARSYSTEM_H
