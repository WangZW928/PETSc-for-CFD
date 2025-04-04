#include <heat.h>

int main(int argc,char **argv) {

    Heat heat;

    PetscCall(PetscInitialize(&argc,&argv,NULL,help));

    heat.Initialze();
    
    heat.Solve();

    PetscCall(PetscFinalize());

    return 0;
}