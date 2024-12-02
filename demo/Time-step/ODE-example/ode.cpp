#include <petsc.h>

class TimeStep {
public:
    TimeStep(PetscInt freedom);
    ~TimeStep();

    PetscErrorCode SetTimeAxis();
    PetscErrorCode Initilze();
    PetscErrorCode SetInitialValuesAndSolve();
    PetscErrorCode ExactSolution();
    static PetscErrorCode FormRHSFunction(TS ts, PetscReal t, Vec y, Vec g, void *ctx);

private:
    Vec y, yexact;
    TS ts;
    PetscInt steps;
    PetscReal t0 = 0.0, tf = 20.0, dt = 0.1, err;
    PetscInt freedom;
};

// Constructor for TimeStep
TimeStep::TimeStep(PetscInt freedom0) {
    freedom = freedom0;
}

PetscErrorCode TimeStep::Initilze() {
    PetscCall(VecCreate(PETSC_COMM_WORLD, &y));
    PetscCall(VecSetSizes(y, PETSC_DECIDE, freedom));
    PetscCall(VecSetFromOptions(y));
    PetscCall(VecDuplicate(y, &yexact));

    PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));
    PetscCall(TSSetProblemType(ts, TS_NONLINEAR));
    PetscCall(TSSetRHSFunction(ts, NULL, TimeStep::FormRHSFunction, NULL));
    PetscCall(TSSetType(ts, TSRK));

    return 0;
}

// Destructor for TimeStep
TimeStep::~TimeStep() {
    PetscErrorCode ierr;
    ierr = VecDestroy(&y); if (ierr) PetscError(PETSC_COMM_SELF, __LINE__, "VecDestroy error", __FILE__, ierr, PETSC_ERROR_REPEAT, " ");
    ierr = VecDestroy(&yexact); if (ierr) PetscError(PETSC_COMM_SELF, __LINE__, "VecDestroy error", __FILE__, ierr, PETSC_ERROR_REPEAT, " ");
    ierr = TSDestroy(&ts); if (ierr) PetscError(PETSC_COMM_SELF, __LINE__, "TSDestroy error", __FILE__, ierr, PETSC_ERROR_REPEAT, " ");
}

// Set the time axis for the solver
PetscErrorCode TimeStep::SetTimeAxis() {

    /*
    
    set the boundary of time axis and time step size
    
    */
    PetscCall(TSSetTime(ts, t0));
    PetscCall(TSSetMaxTime(ts, tf));
    PetscCall(TSSetTimeStep(ts, dt));

    /*
    
    makes sure that adaptive TS types respect the final time you just set.

    */

    PetscCall(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP));
    PetscCall(TSSetFromOptions(ts));
    return 0;
}

// Set initial values and solve
PetscErrorCode TimeStep::SetInitialValuesAndSolve() {
    PetscCall(TSGetTime(ts, &t0));

    PetscReal *ay;
    PetscCall(VecGetArray(y, &ay));
    ay[0] = t0 - PetscSinReal(t0);
    ay[1] = 1.0 - PetscCosReal(t0);
    PetscCall(VecRestoreArray(y, &ay));
    
    PetscCall(TSSolve(ts, y));
    return 0;
}

// Exact solution
PetscErrorCode TimeStep::ExactSolution() {
    PetscCall(TSGetStepNumber(ts, &steps));
    PetscCall(TSGetTime(ts, &tf));

    PetscReal *ay;
    PetscCall(VecGetArray(yexact, &ay));
    ay[0] = tf - PetscSinReal(tf);
    ay[1] = 1.0 - PetscCosReal(tf);
    PetscCall(VecRestoreArray(yexact, &ay));


    PetscCall(VecAXPY(y, -1.0, yexact));  // y <- y - yexact
    PetscCall(VecNorm(y, NORM_INFINITY, &err));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "error at tf = %.3f with %d steps:  |y-y_exact|_inf = %g\n", tf, steps, err));

    return 0;
}

// RHS function
PetscErrorCode TimeStep::FormRHSFunction(TS ts, PetscReal t, Vec y, Vec g, void *ctx) {
    const PetscReal *ay;
    PetscReal *ag;
    PetscCall(VecGetArrayRead(y, &ay));
    PetscCall(VecGetArray(g, &ag));
    ag[0] = ay[1];            // g_1(t,y)
    ag[1] = -ay[0] + t;      // g_2(t,y)
    PetscCall(VecRestoreArrayRead(y, &ay));
    PetscCall(VecRestoreArray(g, &ag));
    return 0;
}

// Main function
int main(int argc, char **argv) {
    PetscCall(PetscInitialize(&argc, &argv, NULL, "ODE system solver example"));

    TimeStep ts(2);

    PetscCall(ts.Initilze());

    PetscCall(ts.SetTimeAxis());
    PetscCall(ts.SetInitialValuesAndSolve());
    PetscCall(ts.ExactSolution());

    PetscCall(PetscFinalize());
    return 0;
}
