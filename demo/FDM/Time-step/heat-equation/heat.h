#include <petsc.h>

class HeatCtx {
public:
    PetscReal D0;    // conductivity
};

class Heat{
public:
    Heat();
    ~Heat();

    HeatCtx user;
    DM da;
    Vec u;
    TS ts;
    DMDALocalInfo info;
    PetscReal t0, tf;
    PetscBool monitorenergy;


    PetscReal f_source(PetscReal x, PetscReal y);
    PetscReal gamma_neumann(PetscReal y);
    PetscErrorCode FormRHSFunctionLocal();
    PetscErrorCode FormRHSJacobianLocal();
    PetscErrorCode Initialze();
    PetscErrorCode Solve();

};

