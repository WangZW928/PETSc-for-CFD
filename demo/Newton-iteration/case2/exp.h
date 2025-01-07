#include <petsc.h>

class Exp {

    public:
        SNES snes;          // nonlinear solver
        Vec x, r;          // solution, residual vectors

        PetscErrorCode Initialize();
        PetscErrorCode Solve();
        PetscErrorCode Destory();
        
        static PetscErrorCode FormFunction(SNES, Vec, Vec, void*);

};