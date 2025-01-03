#include <HeaderFile.h>

PetscErrorCode User :: FormRHSFunctionLocal(Field **aG){

    PetscInt  i,j;
    PetscReal uv2;

    ctx->RHSFcn_called = PETSC_TRUE;

    for (j = info.ys; j < info.ys + info.ym; j++) {
        for (i = info.xs; i < info.xs + info.xm; i++) {
            uv2 = aY[j][i].u * aY[j][i].v * aY[j][i].v;
            aG[j][i].u = - uv2 + ctx->phi * (1.0 - aY[j][i].u);
            aG[j][i].v = + uv2 - (ctx->phi + ctx->kappa) * aY[j][i].v;
        }
    }

}





