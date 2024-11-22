static char help[] = "A structured-grid Poisson solver using DMDA+KSP.\n\n";

#include <petsc.h>

extern PetscErrorCode formMatrix(DM, Mat);
extern PetscErrorCode formExact(DM, Vec);
extern PetscErrorCode formRHS(DM, Vec);

//STARTMAIN
int main(int argc,char **args) {
    DM            da;
    Mat           A;
    Vec           b,u,uexact;
    KSP           ksp;
    PetscReal     errnorm;
    DMDALocalInfo info;

    PetscCall(PetscInitialize(&argc,&args,NULL,help));

    // change default 9x9 size using -da_grid_x M -da_grid_y N
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD,
                 DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,
                 9,9,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&da));
    /*
    PETSC_DECIDE: use default value, means PETSc will decide the value
    dof: degree of freedom, 1 means scalar, 2 means vector, 3 means tensor
    lx,ly: lx[i] means the size of i processor in x direction,null means PETSc will decide
    DM_BOUNDARY_NONE: the value on the boundary is known or can be set, PETSC don't need extra memory to store beyond the boundary points
    DMDA_STENCIL_STAR: the stencil type, STAR means 5-point stencil, BOX means 9-point stencil, the points are (i,j), (i-1,j), (i+1,j), (i,j-1), (i,j+1) used compute the value at (i,j)
    s: the order of points in STENCIL  
    */

    // create linear system matrix A
    PetscCall(DMSetFromOptions(da));
    PetscCall(DMSetUp(da));
    PetscCall(DMCreateMatrix(da,&A));
    PetscCall(MatSetFromOptions(A));

    // create RHS b, approx solution u, exact solution uexact
    PetscCall(DMCreateGlobalVector(da,&b));
    PetscCall(VecDuplicate(b,&u));
    PetscCall(VecDuplicate(b,&uexact));

    // fill vectors and assemble linear system
    PetscCall(formExact(da,uexact));
    PetscCall(formRHS(da,b));
    PetscCall(formMatrix(da,A));

    // create and solve the linear system
    PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));
    PetscCall(KSPSetOperators(ksp,A,A));
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(KSPSolve(ksp,b,u));

    // report on grid and numerical error
    PetscCall(VecAXPY(u,-1.0,uexact));    // u <- u + (-1.0) uxact
    PetscCall(VecNorm(u,NORM_INFINITY,&errnorm));
    PetscCall(DMDAGetLocalInfo(da,&info));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                "on %d x %d grid:  error |u-uexact|_inf = %g\n",
                info.mx,info.my,errnorm));

    PetscCall(VecDestroy(&u));
    PetscCall(VecDestroy(&uexact));
    PetscCall(VecDestroy(&b));
    PetscCall(MatDestroy(&A));
    PetscCall(KSPDestroy(&ksp));
    PetscCall(DMDestroy(&da));
    PetscCall(PetscFinalize());
    return 0;
}
//ENDMAIN

//STARTMATRIX
PetscErrorCode formMatrix(DM da, Mat A) {
    DMDALocalInfo  info;

    /*
    
    typedef struct {
        PetscInt k, j, i, c;
    } MatStencil;

    The i,j, and k represent the logical coordinates over the entire grid 
    (for 2 and 1 dimensional problems the k and j entries are ignored). 
    The c represents the degrees of freedom at each grid point (the dof argument to DMDASetDOF()). 
    If dof is 1 then this entry is ignored.
    
    */



    MatStencil     row, col[5];
    PetscReal      hx, hy, v[5];
    PetscInt       i, j, ncols;

    PetscCall(DMDAGetLocalInfo(da,&info));
    hx = 1.0/(info.mx-1);  hy = 1.0/(info.my-1);
    for (j = info.ys; j < info.ys+info.ym; j++) {
        for (i = info.xs; i < info.xs+info.xm; i++) {
            row.j = j;           // row of A corresponding to (x_i,y_j)
            row.i = i;
            col[0].j = j;        // diagonal entry
            col[0].i = i;
            ncols = 1;
            if (i==0 || i==info.mx-1 || j==0 || j==info.my-1) {
                v[0] = 1.0;      // on boundary: trivial equation
            } else {
                v[0] = 2*(hy/hx + hx/hy); // interior: build a row
                if (i-1 > 0) {
                    col[ncols].j = j;    col[ncols].i = i-1;
                    v[ncols++] = -hy/hx;
                }
                if (i+1 < info.mx-1) {
                    col[ncols].j = j;    col[ncols].i = i+1;
                    v[ncols++] = -hy/hx; // after done this; ncols = nclos + 1
                }
                if (j-1 > 0) {
                    col[ncols].j = j-1;  col[ncols].i = i;
                    v[ncols++] = -hx/hy;
                }
                if (j+1 < info.my-1) {
                    col[ncols].j = j+1;  col[ncols].i = i;
                    v[ncols++] = -hx/hy;
                }
            }
            PetscCall(MatSetValuesStencil(A,1,&row,ncols,col,v,INSERT_VALUES));
        }
    }
    PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
    return 0;
}
//ENDMATRIX

//STARTEXACT
PetscErrorCode formExact(DM da, Vec uexact) {
    PetscInt       i, j;
    PetscReal      hx, hy, x, y, **auexact;
    DMDALocalInfo  info;

    PetscCall(DMDAGetLocalInfo(da,&info));// get local grid info
    hx = 1.0/(info.mx-1);  hy = 1.0/(info.my-1);// grid spacing , info.mx = global x grid size
    PetscCall(DMDAVecGetArray(da, uexact, &auexact));// get array of uexact, 2D array
    for (j = info.ys; j < info.ys+info.ym; j++) {// loop over local grid,info.ys = local y start, info.ym = local y size
        y = j * hy;
        for (i = info.xs; i < info.xs+info.xm; i++) {
            x = i * hx;
            auexact[j][i] = x*x * (1.0 - x*x) * y*y * (y*y - 1.0);// exact solution
        }
    }
    PetscCall(DMDAVecRestoreArray(da, uexact, &auexact));// restore 2D array to Petsc Vec
    return 0;
}
// use Vec is good for parallel computing, but not for 2D array

PetscErrorCode formRHS(DM da, Vec b) {
    PetscInt       i, j;
    PetscReal      hx, hy, x, y, f, **ab;
    DMDALocalInfo  info;

    PetscCall(DMDAGetLocalInfo(da,&info));
    hx = 1.0/(info.mx-1);  hy = 1.0/(info.my-1);
    PetscCall(DMDAVecGetArray(da, b, &ab));
    for (j=info.ys; j<info.ys+info.ym; j++) {
        y = j * hy;
        for (i=info.xs; i<info.xs+info.xm; i++) {
            x = i * hx;
            if (i==0 || i==info.mx-1 || j==0 || j==info.my-1) {
                ab[j][i] = 0.0;  // on boundary: 1*u = 0
            } else {
                f = 2.0 * ( (1.0 - 6.0*x*x) * y*y * (1.0 - y*y)
                    + (1.0 - 6.0*y*y) * x*x * (1.0 - x*x) );
                ab[j][i] = hx * hy * f;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(da, b, &ab));
    return 0;
}
//ENDEXACT
