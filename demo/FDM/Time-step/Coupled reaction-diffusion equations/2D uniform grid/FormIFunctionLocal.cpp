#include <HeaderFile.h>

PetscErrorCode User::FormIFunctionLocal(Field **aYdot, Field **aF){

    PetscInt         i, j;
    const PetscReal  h = ctx.L / (PetscReal)(info.mx),
                    Cu = ctx.Du / (6.0 * h * h),
                    Cv = ctx.Dv / (6.0 * h * h);
    PetscReal       u, v, lapu, lapv;

    ctx.IFcn_called = PETSC_TRUE;

    for (j = info->ys; j < info->ys + info->ym; j++) {
        for (i = info->xs; i < info->xs + info->xm; i++) {
            u = aY[j][i].u;
            v = aY[j][i].v;
            lapu =     aY[j+1][i-1].u + 4.0*aY[j+1][i].u +   aY[j+1][i+1].u
                 + 4.0*aY[j][i-1].u -    20.0*u        + 4.0*aY[j][i+1].u
                 +   aY[j-1][i-1].u + 4.0*aY[j-1][i].u +   aY[j-1][i+1].u;
            lapv =     aY[j+1][i-1].v + 4.0*aY[j+1][i].v +   aY[j+1][i+1].v
                 + 4.0*aY[j][i-1].v -    20.0*v        + 4.0*aY[j][i+1].v
                 +   aY[j-1][i-1].v + 4.0*aY[j-1][i].v +   aY[j-1][i+1].v;
            aF[j][i].u = aYdot[j][i].u - Cu * lapu;
            aF[j][i].v = aYdot[j][i].v - Cv * lapv;
        }
    }
  return 0;

}