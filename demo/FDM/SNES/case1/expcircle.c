//STARTWHOLE
static char help[] = "Newton's method for a two-variable system.\n"
"No analytical Jacobian.  Run with -snes_fd or -snes_mf.\n\n";

#include <petsc.h>

extern PetscErrorCode FormFunction(SNES, Vec, Vec, void*);

int main(int argc,char **argv) {
    SNES  snes;          // nonlinear solver
    Vec   x, r;          // solution, residual vectors

    PetscCall(PetscInitialize(&argc,&argv,NULL,help));
    PetscCall(VecCreate(PETSC_COMM_WORLD,&x));
    PetscCall(VecSetSizes(x,PETSC_DECIDE,2)); // set size of x, where PETSC_DECIDE means local size is determined by the solver
    PetscCall(VecSetFromOptions(x));
    PetscCall(VecSet(x,1.0));         // initial iterate
    PetscCall(VecDuplicate(x,&r));

    PetscCall(SNESCreate(PETSC_COMM_WORLD,&snes));//create a nonlinear solver
    PetscCall(SNESSetFunction(snes,r,FormFunction,NULL));//set the function evaluation routine for error calculation
    PetscCall(SNESSetFromOptions(snes));//allow command line options to be passed to the solver and how the Jacobian is calculated
    /*
    
    -snes_fd: use finite difference to calculate the Jacobian
    -snes_mf: not assembled, but, in a Krylov iterative method for solving system , the action of the Jacobian on vectors (i.e., JF(xk)y) is computed by finite differences
    -snes_fd_color: computed and assembled by calling FormFunction() substantially fewer than N times to compute F(xk + δv) for special vectors v, by using a graph-coloring algorithm based on the Jacobian sparsity pattern to construct the vectors v

    */
    PetscCall(SNESSolve(snes,NULL,x));// store the solution in x, NULL means no rhs
    PetscCall(VecView(x,PETSC_VIEWER_STDOUT_WORLD));

    PetscCall(SNESDestroy(&snes));
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&r));
    PetscCall(PetscFinalize());
    return 0;
}

PetscErrorCode FormFunction(SNES snes, Vec x, Vec F, void *ctx) {
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
//ENDWHOLE
