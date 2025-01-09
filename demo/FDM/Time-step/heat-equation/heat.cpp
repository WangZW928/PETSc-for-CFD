#include <heat.h>

Heat::Heat() {

    user.D0 = 1.0;
    monitorenergy = PETSC_FALSE;

}

Heat::~Heat() {

    PetscCall(TSDestroy(&ts));
    PetscCall(DMDestroy(&da));
    PetscCall(VecDestroy(&u));

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
           "heat.c",user.D0,&user.D0,NULL));
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
    PetscCall(DMDAGetLocalInfo(da,&info));
    this->GetSpaces();

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
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
            "Heat equation on a 2D grid, %D x %D\n"
            "Diffusivity D0 = %g\n"
            "t0 = %g, tf = %g\n",
            info.mx,info.my,(double)user.D0,(double)t0,(double)tf));

    return 0;


}

PetscErrorCode Heat::Solve(){

    PetscCall(VecSet(u,0.0));
    PetscCall(TSSolve(ts,u));

    return 0;

}

PetscErrorCode Heat ::GetSpaces(){

    hx = 1.0/(PetscReal)(info.mx-1);
    hy = 1.0/(PetscReal)(info.my);

    return 0;


}

PetscErrorCode Heat ::FormRHSFunctionLocal(PetscReal t, PetscReal **au, PetscReal **aG){

    PetscInt  i,j,mx = info.mx;
    PetscReal x,y,ul,ur,uxx,uyy;
    PetscInt xs,ys,xe,ye;
    xs = info.xs, ys = info.ys, xe = info.xs + info.xm, ye = info.ys + info.ym;

    for(j = ys;j<ye;j++){
        y = hy*j;
        for(i=xs;i<xe;i++){
            x = hx*i;
            ul = (i == 0) ? au[j][i+1] + 2.0*hx*gamma_neumann(y) : au[j][i-1];//which is used for calcute uxx[j][0]
            ur = (i == mx-1) ? au[j][i-1] : au[j][i+1];//which is used for calcute uxx[j][mx-1],where gamma is 0
            uxx = (ul - 2.0*au[j][i] + ur)/(hx*hx);
            uyy = (au[j-1][i] - 2.0*au[j][i] + au[j+1][i])/(hy*hy);//y is periodic
            aG[j][i] = user.D0*(uxx + uyy) + f_source(x,y);
        }
    }

    return 0;

}

PetscErrorCode Heat ::FormRHSJacobianLocal(PetscReal t, PetscReal **au,Mat J, Mat P){

    // calculate the Jacobian matrix for G(u) = 0
    
    PetscInt i,j,ncols;
    const PetscReal D = user.D0;
    PetscInt max = info.mx;
    PetscInt xs,ys,xe,ye;
    xs = info.xs, ys = info.ys, xe = info.xs + info.xm, ye = info.ys + info.ym;
    PetscReal        hx, hy, hx2, hy2, v[5];
    MatStencil       col[5],row;

    hx2 = hx * hx;  hy2 = hy * hy;

    for(j=ys;j<ye;j++){
        row.j = j;  col[0].j = j;
        for(i=xs;i<xe;i++){
            row.i = i;
            col[0].i = i;
            v[0] = - 2.0 * D * (1.0 / hx2 + 1.0 / hy2);

            col[1].j = j-1;  col[1].i = i;    v[1] = D / hy2;
            col[2].j = j+1;  col[2].i = i;    v[2] = D / hy2;
            col[3].j = j;    col[3].i = i-1;  v[3] = D / hx2;
            col[4].j = j;    col[4].i = i+1;  v[4] = D / hx2;
            ncols = 5;
            if(i==0){
                ncols = 4;
                col[3].j = j;    col[3].i = i+1;  v[3] = D / hx2;
            }else if(i==max-1){
                ncols = 4;
                col[4].j = j;    col[4].i = i-1;  v[4] = D / hx2;
            }
            PetscCall(MatSetValuesStencil(P,1,&row,ncols,col,v,INSERT_VALUES));
        }
    }
    PetscCall(MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY));

    if(J!=P){
        PetscCall(MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY));
    }

    return 0;
}

