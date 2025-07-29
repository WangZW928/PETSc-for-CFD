#include <stdio.h>
#include <stdlib.h>
#include "HYPRE.h"
#include "HYPRE_utilities.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"
#include <mpi.h>


int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);     // ✅ 初始化 MPI

    HYPRE_IJMatrix A;
    HYPRE_ParCSRMatrix par_A;
    int ilower = 0, iupper = 2;
    int jlower = 0, jupper = 2;
    int rows = iupper - ilower + 1;

    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, jlower, jupper, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A);

    // 创建一个简单 3x3 矩阵：
    // [ 2 -1  0 ]
    // [-1  2 -1 ]
    // [ 0 -1  2 ]
    for (int i = ilower; i <= iupper; i++) {
        int ncols = 0;
        int cols[3];
        double values[3];

        if (i > ilower) {
            cols[ncols] = i - 1;
            values[ncols++] = -1.0;
        }

        cols[ncols] = i;
        values[ncols++] = 2.0;

        if (i < iupper) {
            cols[ncols] = i + 1;
            values[ncols++] = -1.0;
        }

        HYPRE_IJMatrixSetValues(A, 1, &ncols, &i, cols, values);
    }

    HYPRE_IJMatrixAssemble(A);
    HYPRE_IJMatrixGetObject(A, (void**) &par_A);

    // 打印矩阵（行主序）
    printf("==> CSR 格式输出:\n");
    for (int i = ilower; i <= iupper; i++) {
        int ncols = 3;
        int cols[3];
        double values[3];

        HYPRE_IJMatrixGetValues(A, 1, &ncols, &i, cols, values);

        printf("row %d: ", i);
        for (int j = 0; j < ncols; j++) {
            printf("(%d, %.1f) ", cols[j], values[j]);
        }
        printf("\n");
    }

    HYPRE_IJMatrixDestroy(A);
    MPI_Finalize();               // ✅ 结束 MPI
    return 0;
}

