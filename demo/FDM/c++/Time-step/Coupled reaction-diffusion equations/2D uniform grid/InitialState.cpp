#include <HeaderFile.h>

PetscErrorCode User::InitialState() {
    
    PetscInt         i,j;
    PetscReal        sx,sy;
    const PetscReal  ledge = (ctx.L - 0.5) / 2.0, // nontrivial initial values on
                    redge = ctx.L - ledge;       //   ledge < x,y < redge
    DMDACoor2d       **aC;
    Field            **aY;

    PetscCall(VecSet(this.x,0.0));

    if (this.noiselevel > 0.0) {
      // noise added to usual initial condition is uniform on [0,noiselevel],
      //     independently for each location and component
        PetscCall(VecSetRandom(this.x,NULL));
        PetscCall(VecScale(this.x,this.noiselevel));
      
    }

    PetscCall(DMDAGetLocalInfo(this.da,&this.info));
    PetscCall(DMDAGetCoordinateArray(this.da,&aC));
    PetscCall(DMDAVecGetArray(this.da,this.x,&aY));

    for (j = this.info.ys; j < this.info.ys+this.info.ym; j++) {
        for (i = this.info.xs; i < this.info.xs+this.info.xm; i++) {

            if ((aC[j][i].x >= ledge) && (aC[j][i].x <= redge)
                && (aC[j][i].y >= ledge) && (aC[j][i].y <= redge)) {

                    sx = PetscSinReal(4.0 * PETSC_PI * aC[j][i].x);
                    sy = PetscSinReal(4.0 * PETSC_PI * aC[j][i].y);
                    aY[j][i].v += 0.5 * sx * sx * sy * sy;
            }
            aY[j][i].u += 1.0 - 2.0 * aY[j][i].v;
        }
    }
  PetscCall(DMDAVecRestoreArray(this.da,this.x,&aY));
  PetscCall(DMDARestoreCoordinateArray(this.da,&aC));
  return 0;

  }

