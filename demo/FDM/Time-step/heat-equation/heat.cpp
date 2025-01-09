#include <heat.h>

Heat::Heat() {

    user.D0 = 1.0;
    monitorenergy = PETSC_FALSE;

}

PetscErrorCode Heat::f_source(PetscReal x, PetscReal y) {
    return 3.0 * PetscExpReal(-25.0 * (x-0.6) * (x-0.6))
               * PetscSinReal(2.0*PETSC_PI*y);
}

PetscErrorCode Heat::gamma_neumann(PetscReal y) {
    return PetscSinReal(6.0 * PETSC_PI * y);
}

PetscErrorCode Heat::Initialze(){

    PetscOptionsBegin(PETSC_COMM_WORLD, "ht_", "options for heat", "");
    PetscCall(PetscOptionsReal("-D0","constant thermal diffusivity",
           "heat.c",D0,&D0,NULL));
    PetscCall(PetscOptionsBool("-monitor","also display total heat energy at each step",
              "heat.c",monitorenergy,&monitorenergy,NULL));
    PetscOptionsEnd();

    PetscCall(
        DMDACreate2d(
            PETSC_COMM_WORLD,
            DM_BOUNDARY_NONE,
            DM_BOUNDARY_PERIODIC,
            DMDA_STENCIL_STAR,
            5,
            4,
            PETSC_DECIDE,
            PETSC_DECIDE,
            1,
            1,
            NULL,
            NULL,
            &da
        )
    );

    PetscCall(DMSetFromOptions(da));//set da from options by user
    PetscCall(DMSetUp(da));//always call after DMSetFromOptions
    PetscCall(DMCreateGlobalVector(da,&u));//create global vector

    PetscCall(TSCreate(PETSC_COMM_WORLD,&ts));
    PetscCall(TSSetProblemType(ts,TS_NONLINEAR));//OR TS_LINEAR
    PetscCall(TSSetDM(ts,da));

    PetscCall(TSSetApplicationContext(ts,&user));//set user context for ts,if we need to pass any data to ts

    PetscCall(DMDATSSetRHSFunctionLocal(da,INSERT_VALUES,
           (DMDATSRHSFunctionLocal)FormRHSFunctionLocal,&user));
    PetscCall(DMDATSSetRHSJacobianLocal(da,
           (DMDATSRHSJacobianLocal)FormRHSJacobianLocal,&user));


    PetscCall(TSSetType(ts,TSBDF));
    PetscCall(TSSetTime(ts,0.0));//Allows one to reset the time.
    PetscCall(TSSetMaxTime(ts,0.1));//Sets the maximum (or final) time for timestepping.
    PetscCall(TSSetTimeStep(ts,0.001));
    PetscCall(TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP));
    PetscCall(TSSetFromOptions(ts));

    PetscCall(TSGetTime(ts,&t0));
    PetscCall(TSGetMaxTime(ts,&tf));
    PetscCall(DMDAGetLocalInfo(da,&info));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
            "Heat equation on a 2D grid, %D x %D\n"
            "Diffusivity D0 = %g\n"
            "t0 = %g, tf = %g\n",
            info.mx,info.my,(double)user.D0,(double)t0,(double)tf));
            



}
