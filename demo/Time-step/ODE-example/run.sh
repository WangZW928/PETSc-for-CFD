#!/bin/bash

# Initialize PETSc environment variables
export PETSC_DIR=/home/wangzewei/PETSc/petsc
export PETSC_ARCH=arch-linux2-c-debug
export PATH=$PETSC_DIR/$PETSC_ARCH/bin:$PATH
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
export PETSC_LIB_DIR=$PETSC_DIR/$PETSC_ARCH/lib

# Read parameters from the file, ignoring comments and empty lines
params=$(grep -v '^\s*#' params.txt | grep -v '^\s*$')

# Run the ODE program with the parameters
./ode $params
