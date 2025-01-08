#ifndef included_Grid
#define included_Grid

#include <petsc.h>
#include <string>
#include <vector>

class Grid {
public:

Grid();
~Grid();

PetscInt Dim;//1,2,3
PetscInt Is1D;//1 or 0
PetscInt Is2D;//1 or 0
PetscInt Is3D;

PetscInt Nx_nodes, Ny_nodes, Nz_nodes; //number of nodes in each direction


PetscErrorCode ReadGrid();
PetscErrorCode ReadGrid1D();
PetscErrorCode ReadGrid2D();
PetscErrorCode ReadGrid3D();

PetscErrorCode CreateGrid();
PetscErrorCode CreateGrid1D();
PetscErrorCode CreateGrid2D();
PetscErrorCode CreateGrid3D();

private:


};





#endif



