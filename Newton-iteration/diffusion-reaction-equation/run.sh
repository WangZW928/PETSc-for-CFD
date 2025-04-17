# Read parameters from the file, ignoring comments and empty lines
params=$(grep -v '^\s*#' params.txt | grep -v '^\s*$')

# Run the ODE program with the parameters
mpirun -n 1 ./test $params