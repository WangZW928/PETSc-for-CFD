#include "exp.h"  // 使用双引号而不是尖括号


PetscErrorCode Exp :: Initialize(){
    PetscCall(VecCreate(PETSC_COMM_WORLD, &x));
    PetscCall(VecSetSizes(x, PETSC_DECIDE, 2));//set size of x, where PETSC_DECIDE means local size is determined by the solver
    PetscCall(VecSetFromOptions(x));
    PetscCall(VecSet(x, 1.0));         // initial iterate
    PetscCall(VecDuplicate(x, &r));

    PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));//create a nonlinear solver
    PetscCall(SNESSetFromOptions(snes));//allow command line options to be passed to the solver and how the Jacobian is calculated

    return 0;
}

PetscErrorCode Exp::Solve(){
    PetscCall(SNESSetFunction(snes, r, Exp::FormFunction, NULL));//set the function evaluation routine for error calculation
    PetscCall(SNESSolve(snes, NULL, x));// store the solution in x, NULL means no rhs
    PetscCall(VecView(x,PETSC_VIEWER_STDOUT_WORLD));

    return 0;

}

PetscErrorCode Exp::Destory(){
    PetscCall(SNESDestroy(&snes));
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&r));
    return 0;
}

PetscErrorCode Exp::FormFunction(SNES snes, Vec x, Vec F, void *ctx){
    const PetscReal  b = 2.0, *ax;
    PetscReal        *aF;

    PetscCall(VecGetArrayRead(x,&ax));

    PetscCall(VecGetArray(F,&aF));

    aF[0] = (1.0 / b) * PetscExpReal(b * ax[0]) - ax[1];
    aF[1] = ax[0] * ax[0] + ax[1] * ax[1] - 1.0;
    PetscCall(VecRestoreArrayRead(x,&ax));
    PetscCall(VecRestoreArray(F,&aF));
    return 0;

}