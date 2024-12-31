#include <petsc.h>

static char help[] =
"Coupled reaction-diffusion equations (Pearson 1993).  Option prefix -ptn_.\n"
"Demonstrates form  F(t,Y,dot Y) = G(t,Y)  where F() is IFunction and G() is\n"
"RHSFunction().  Implements IJacobian() and RHSJacobian().  Defaults to\n"
"ARKIMEX (= adaptive Runge-Kutta implicit-explicit) TS type.\n\n";

class Field{

  public:
    PetscReal u, v;

};

class PatternCtx{
  public:
    PatternCtx();

    PetscReal  L,     // domain side length
             Du,    // diffusion coefficient: u equation
             Dv,    //                        v equation
             phi,   // "dimensionless feed rate" (F in Pearson 1993)
             kappa; // "dimensionless rate constant" (k in Pearson 1993)
    PetscBool  IFcn_called, IJac_called, RHSFcn_called, RHSJac_called;

    void PetscOption();

};


class User{
public:
  TS        ts;
  Vec       x;
  DM        da;

  DMDALocalInfo  info;
  PetscReal      noiselevel = -1.0;  // negative value means no initial noise
  PetscBool      no_rhsjacobian = PETSC_FALSE,
                 no_ijacobian = PETSC_FALSE,
                 call_back_report = PETSC_FALSE;
  TSType         type;

  PatternCtx     ctx;

  void PetscOption();
  PetscErrorCode Initilze();
  PetscErrorCode STARTTSSETUP();


  PetscErrorCode InitialState();
  PetscErrorCode FormRHSFunctionLocal();
  PetscErrorCode FormRHSJacobianLocal();
  PetscErrorCode FormIFunctionLocal();
  PetscErrorCode FormIJacobianLocal();
};
