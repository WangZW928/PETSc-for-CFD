#include <HeaderFile.h>

PatternCtx :: PatternCtx(){
    L = 2.5;
    Du = 8.0e-5;
    Dv = 4.0e-5;
    phi = 0.024;
    kappa = 0.06;
    IFcn_called = PETSC_FALSE;
    IJac_called = PETSC_FALSE;
    RHSFcn_called = PETSC_FALSE;
    RHSJac_called = PETSC_FALSE;
}

void PatternCtx::PetscOption(){

  PetscOptionsBegin(PETSC_COMM_WORLD, "ptn_", "options for patterns", "");

  PetscCall(PetscOptionsReal("-Du","diffusion coefficient of first equation",
           "pattern.c",Du,&Du,NULL));

  PetscCall(PetscOptionsReal("-Dv","diffusion coefficient of second equation",
           "pattern.c",Dv,&Dv,NULL));

  PetscCall(PetscOptionsReal("-kappa","dimensionless rate constant (=k in (Pearson, 1993))",
           "pattern.c",kappa,&kappa,NULL));

  PetscCall(PetscOptionsReal("-L","square domain side length; recommend L >= 0.5",
           "pattern.c",L,&L,NULL));

  PetscCall(PetscOptionsReal("-phi","dimensionless feed rate (=F in (Pearson, 1993))",
           "pattern.c",phi,&phi,NULL));
  PetscOptionsEnd();

  return 0;
}

void User::PetscOption(){

  PetscOptionsBegin(PETSC_COMM_WORLD, "ptn_", "options for patterns", "");

  PetscCall(PetscOptionsBool("-call_back_report","report on which user-supplied call-backs were actually called",
           "pattern.c",call_back_report,&(call_back_report),NULL));

  
  PetscCall(PetscOptionsBool("-no_ijacobian","do not set call-back DMDATSSetIJacobian()",
           "pattern.c",no_ijacobian,&(no_ijacobian),NULL));

  PetscCall(PetscOptionsBool("-no_rhsjacobian","do not set call-back DMDATSSetRHSJacobian()",
           "pattern.c",no_rhsjacobian,&(no_rhsjacobian),NULL));

  PetscCall(PetscOptionsReal("-noisy_init",
           "initialize u,v with this much random noise (e.g. 0.2) on top of usual initial values",
           "pattern.c",noiselevel,&noiselevel,NULL));

  PetscOptionsEnd();

  return 0;
}

PetscErrorCode User::Initilze(){

    ctx.PetscOption();

    this->PetscOption();

    PetscCall(DMDACreate2d(PETSC_COMM_WORLD,
               DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
               DMDA_STENCIL_BOX,  // for 9-point stencil
               3,3,PETSC_DECIDE,PETSC_DECIDE,
               2, 1,              // degrees of freedom, stencil width
               NULL,NULL,&da));

    PetscCall(DMSetFromOptions(da));
    PetscCall(DMSetUp(da));
    PetscCall(DMDASetFieldName(da,0,"u"));
    PetscCall(DMDASetFieldName(da,1,"v"));
    PetscCall(DMDAGetLocalInfo(da,&info));


    if (info.mx != info.my) {
      SETERRQ(PETSC_COMM_SELF,1,"code  requires mx == my");
    }


    PetscCall(DMDASetUniformCoordinates(da, 0.0, ctx.L, 0.0, ctx.L, -1.0, -1.0));


    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
           "running on %d x %d grid with square cells of side h = %.6f ...\n",
           info.mx,info.my,ctx.L/(PetscReal)(info.mx)));

    



}

PetscErrorCode User::STARTTSSETUP(){

    PetscCall(TSCreate(PETSC_COMM_WORLD,&ts));
    PetscCall(TSSetProblemType(ts,TS_NONLINEAR));
    PetscCall(TSSetDM(ts,da));
    PetscCall(TSSetApplicationContext(ts,&ctx));

    PetscCall(DMDATSSetRHSFunctionLocal(da,INSERT_VALUES,
           (DMDATSRHSFunctionLocal)this.FormRHSFunctionLocal,&ctx));

    if (!no_rhsjacobian) {
      PetscCall(DMDATSSetRHSJacobianLocal(da,
               (DMDATSRHSJacobianLocal)this.FormRHSJacobianLocal,&ctx));
    }

    PetscCall(DMDATSSetIFunctionLocal(da,INSERT_VALUES,
           (DMDATSIFunctionLocal)this.FormIFunctionLocal,&ctx));
    if (!no_ijacobian) {
        PetscCall(DMDATSSetIJacobianLocal(da,
               (DMDATSIJacobianLocal)this.FormIJacobianLocal,&ctx));
    }

    PetscCall(TSSetType(ts,TSARKIMEX));
    PetscCall(TSSetTime(ts,0.0));
    PetscCall(TSSetMaxTime(ts,200.0));
    PetscCall(TSSetTimeStep(ts,5.0));
    PetscCall(TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP));
    PetscCall(TSSetFromOptions(ts));

}

