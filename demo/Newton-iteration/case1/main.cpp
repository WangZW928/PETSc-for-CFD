#include "exp.h"  // 使用双引号而不是尖括号


int main(int argc, char **argv) {
    PetscCall(PetscInitialize(&argc, &argv, NULL, "ODE system solver example"));

    Exp exp;
    PetscCall(exp.Initialize());
    PetscCall(exp.Solve());
    PetscCall(exp.Destory());

    PetscCall(PetscFinalize());
    
    return 0;
}


    /*
    -snes_fd: use finite difference to calculate the Jacobian
    -snes_mf: not assembled, but, in a Krylov iterative method for solving system , the action of the Jacobian on vectors (i.e., JF(xk)y) is computed by finite differences
    -snes_fd_color: computed and assembled by calling FormFunction() substantially fewer than N times to compute F(xk + δv) for special vectors v, by using a graph-coloring algorithm based on the Jacobian sparsity pattern to construct the vectors v

    */