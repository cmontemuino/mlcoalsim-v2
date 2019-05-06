#!/bin/bash -
#title          :compile.sh
#description    :This script compiles and generates all of the ml_coalsim binaries

mpicc source/*.c -lm -o bin/mlcoalsimXmpi -Wall -pedantic -DBENCHMARK -DinMPI -O3
mpicc source/*.c -lm -o bin/mlcoalsimXmpi_ZnS -Wall -pedantic -DBENCHMARK -DinMPI -DZNS_ACTIVE -O3
gcc source/*.c -lm -o bin/mlcoalsimX_ZnS -Wall -pedantic -DBENCHMARK -DZNS_ACTIVE -O3
gcc source/*.c -lm -o bin/mlcoalsimX -Wall -pedantic -DBENCHMARK -O3
