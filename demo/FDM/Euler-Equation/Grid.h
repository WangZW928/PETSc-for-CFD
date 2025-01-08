#ifndef included_Grid
#define included_Grid

#include <petsc.h>
#include <string>
#include <vector>

class Field{
    public:
    PetscScalar x,y,z;
}


class Grid {
public:

Grid();
~Grid();

DM da_s;//scatter context
DM da_v;//vector context
PetscInt Dim;//1,2,3
PetscInt Is1D;//1 or 0
PetscInt Is2D;//1 or 0
PetscInt Is3D;
PetscInt dof;//number of degrees of freedom per node
PetscInt s;//number of stencil points

PetscInt i_periodic;//1 or 0
PetscInt j_periodic;
PetscInt k_periodic;
//if periodic, then the boundary is periodic

PetscInt Nx_nodes, Ny_nodes, Nz_nodes; //number of nodes in each direction


PetscErrorCode ReadGrid();
PetscErrorCode ReadGrid1D();
PetscErrorCode ReadGrid2D();
PetscErrorCode ReadGrid3D();

PetscErrorCode CreateGrid1D();
PetscErrorCode CreateGrid2D();
PetscErrorCode CreateGrid3D();

private:


};







#endif



