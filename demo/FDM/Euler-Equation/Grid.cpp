#include <Grid.h>

PetscErrorCode Grid::CreateGrid1D(){

}

PetscErrorCode Grid::CreateGrid2D(){

}

PetscErrorCode Grid::CreateGrid3D(){

    PetscInt size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    PetscInt m, n, p;
    DMDABoundaryType bx=DMDA_BOUNDARY_GHOSTED, 
                     by=DMDA_BOUNDARY_GHOSTED, 
                     bz=DMDA_BOUNDARY_GHOSTED;
    m = n = p = PETSC_DECIDE;
        
    if (i_periodic) bx = DMDA_BOUNDARY_PERIODIC;
    if (j_periodic) by = DMDA_BOUNDARY_PERIODIC;
    if (k_periodic) bz = DMDA_BOUNDARY_PERIODIC;

    DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_BOX,
                 d_IM+1, d_JM+1, d_KM+1, m, n,
                 p, 1, s, PETSC_NULL, PETSC_NULL, PETSC_NULL,
                 &da_s);

    DMDAGetInfo(da_s, PETSC_NULL, PETSC_NULL, PETSC_NULL, 
                PETSC_NULL, &m, &n, &p, PETSC_NULL, PETSC_NULL, 
                PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);
    PetscPrintf(PETSC_COMM_WORLD, "**DM 3D Proc Distribution: %i %i %i\n", 
                                  m, n, p);

    DMDASetUniformCoordinates(da_s, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

    DMGetCoordinateDM(da_s, &da_v);

    return 0;




}

PetscErrorCode Grid::ReadGrid(){

    MPI_Barrier(PETSC_COMM_WORLD);
    PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
    PetscOptionsGetReal(PETSC_NULL, "-dim", &Dim, PETSC_NULL);

    switch(Dim){
        case 1:
            Is1D = 1;
            Is2D = 0;
            Is3D = 0;

            PetscPrintf(PETSC_COMM_WORLD, "--- Solve 1D problem ---\n");

            break;
        case 2:
            Is1D = 0;
            Is2D = 1;
            Is3D = 0;

            PetscPrintf(PETSC_COMM_WORLD, "--- Solve 2D problem ---\n");

            break;
        case 3:
            Is1D = 0;
            Is2D = 0;
            Is3D = 1;

            PetscPrintf(PETSC_COMM_WORLD, "--- Solve 3D problem ---\n");
            break;
    }

    
    char str[256];// absolute path to the file
    sprintf(str, "%s", "xyz.dat");

    if(Is1D){
        ReadGrid1D(str);
    }
    else if(Is2D){
        ReadGrid2D(str);
    }
    else if(Is3D){
        ReadGrid3D(str);
    }



}


PetscErrorCode Grid::ReadGrid1D(char &str){

    FILE *fp;
    fd = fopen(str, "r");
    PetscPrintf(PETSC_COMM_WORLD, "Reading: %s\n", str);
    if (fd==NULL) 
    printf("Cannot open %s !\n", str),exit(0);



}

PetscErrorCode Grid::ReadGrid2D(char &str){

    FILE *fp;
    fd = fopen(str, "r");
    PetscPrintf(PETSC_COMM_WORLD, "Reading: %s\n", str);
    if (fd==NULL) 
    printf("Cannot open %s !\n", str),exit(0);



}

PetscErrorCode Grid::ReadGrid3D(char &str){

    FILE *fp;
    fd = fopen(str, "r");
    PetscPrintf(PETSC_COMM_WORLD, "Reading: %s\n", str);
    if (fd==NULL) 
    printf("Cannot open %s !\n", str),exit(0);

    std::vector<double> X, Y, Z;
    double tmp;

    PetscInt i,j,k;

    fscanf(fd, "%i %i %i\n", &Nx_nodes, &Ny_nodes, &Nz_nodes);
    X.resize(Nx_nodes);
    Y.resize(Ny_nodes);
    Z.resize(Nz_nodes);

    PetscPrintf(PETSC_COMM_WORLD, "Reading %s %dx%dx%d\n", 
                                       str, Nx_nodes, Ny_nodes, Nz_nodes);

    for (i=0; i<Nx_nodes; i++) fscanf(fd, "%le %le %le\n", &X[i], &tmp, &tmp);
    for (j=0; j<Ny_nodes; j++) fscanf(fd, "%le %le %le\n", &tmp, &Y[j], &tmp);
    for (k=0; k<Nz_nodes; k++) fscanf(fd, "%le %le %le\n", &tmp, &tmp, &Z[k]);

    CreateGrid3D(); 
    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Created DM\n");

    DMDALocalInfo info;

    /*

typedef struct {
  PetscInt        dim, dof, sw;
  PetscInt        mx, my, mz;     // global number of grid points in each direction 
  PetscInt        xs, ys, zs;     // starting point of this processor, excluding ghosts 
  PetscInt        xm, ym, zm;     // number of grid points on this processor, excluding ghosts 
  PetscInt        gxs, gys, gzs;  // starting point of this processor including ghosts 
  PetscInt        gxm, gym, gzm;  // number of grid points on this processor including ghosts 
  DMBoundaryType  bx, by, bz;     // type of ghost nodes at boundary 
  DMDAStencilType st;
  DM              da;
} DMDALocalInfo;
    
    */

    DMDAGetLocalInfo(da_s, &info);

    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;

    Field ***coord;
    Vec Coord, gCoord;
    DMGetCoordinatesLocal(da_s, &Coord);
    DMDAVecGetArray(da_v, Coord, &coord);


    for (k=zs; k<ze; k++) {
        for (j=ys; j<ye; j++) {
            for (i=xs; i<xe; i++) {

                if(k>=0 && k<Nz_nodes && j>=0 && j<Ny_nodes && i>=0 && i<Nx_nodes){
                    coord[k][j][i].x = X[i];
                    coord[k][j][i].y = Y[j];
                    coord[k][j][i].z = Z[k];
                }

            }
        }
    }

    DMDAVecRestoreArray(da_v, Coord, &coord);
    DMgetCoordinates(da_s, &gCoord);

    DMLocalToGlobalBegin(da_s, Coord, INSERT_VALUES, gCoord);
    DMLocalToGlobalEnd(da_s, Coord, INSERT_VALUES, gCoord);

    DMGlobalToLocalBegin(da_s, gCoord, INSERT_VALUES, Coord);
    DMGlobalToLocalEnd(da_s, gCoord, INSERT_VALUES, Coord);

    fclose(fd);

    MPI_Barrier(PETSC_COMM_WORLD);

    return 0;

}

Grid::Grid(){

}

Grid::~Grid(){

}

