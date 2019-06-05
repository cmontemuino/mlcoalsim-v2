#!/bin/bash -
#title          :compile.sh
#description    :This script compiles and generates all of the ml_coalsim binaries

LOG_LEVEL=ZF_LOG_INFO
BENCHMARK_ENABLED=0

mpicc source/*.c -lm -o bin/mlcoalsimXmpi -Wall -pedantic -DBENCHMARK_ENABLED=$BENCHMARK_ENABLED -DZF_LOG_DEF_LEVEL=$LOG_LEVEL -DinMPI -O3
mpicc source/*.c -lm -o bin/mlcoalsimXmpi_ZnS -Wall -pedantic -DBENCHMARK_ENABLED=$BENCHMARK_ENABLED -DZF_LOG_DEF_LEVEL=$LOG_LEVEL -DinMPI -DZNS_ACTIVE -O3
gcc source/*.c -lm -o bin/mlcoalsimX_ZnS -Wall -pedantic -DBENCHMARK_ENABLED=$BENCHMARK_ENABLED -DZF_LOG_DEF_LEVEL=$LOG_LEVEL -DZNS_ACTIVE -O3
gcc source/*.c -lm -o bin/mlcoalsimX -Wall -pedantic -DBENCHMARK_ENABLED=$BENCHMARK_ENABLED -DZF_LOG_DEF_LEVEL=$LOG_LEVEL -O3
