#include <Grid.h>

PetscErrorCode Grid::CreateGrid(){

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





}

Grid::Grid(){

}

Grid::~Grid(){

}

