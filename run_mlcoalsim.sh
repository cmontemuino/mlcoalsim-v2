echo sh ./run_mlcoalsim.sh
echo cd source
cd source
echo gcc *.c -lm -o mlcoalsimX_ZnS -Wall -pedantic -DZNS_ACTIVE -O3
gcc *.c -lm -o mlcoalsimX_ZnS -Wall -pedantic -DZNS_ACTIVE -O3
echo gcc *.c -lm -o mlcoalsimX -Wall -pedantic -O3
gcc *.c -lm -o mlcoalsimX -Wall -pedantic -O3
echo mpicc *.c -lm -o mlcoalsimXmpi_ZnS -Wall -pedantic -DZNS_ACTIVE -O3 -DinMPI
mpicc *.c -lm -o mlcoalsimXmpi_ZnS -Wall -pedantic -DZNS_ACTIVE -O3 -DinMPI
echo mpicc *.c -lm -o mlcoalsimXmpi -Wall -pedantic -O3 -DinMPI
mpicc *.c -lm -o mlcoalsimXmpi -Wall -pedantic -O3 -DinMPI
echo mv mlcoalsimX_ZnS ../bin
mv mlcoalsimX_ZnS ../bin
echo mv mlcoalsimX ../bin
mv mlcoalsimX ../bin
echo mv mlcoalsimXmpi_ZnS ../bin
mv mlcoalsimXmpi_ZnS ../bin
echo mv mlcoalsimXmpi ../bin
mv mlcoalsimXmpi ../bin
echo cd ../bin
cd ../bin
echo cp ./mlcoalsimX ../examples/example00
cp ./mlcoalsimX ../examples/example00
echo cp ./mlcoalsimX_ZnS ../examples/example00
cp ./mlcoalsimX_ZnS ../examples/example00
echo cp ./mlcoalsimXmpi ../examples/example00
cp ./mlcoalsimXmpi ../examples/example00
echo cp ./mlcoalsimXmpi_ZnS ../examples/example00
cp ./mlcoalsimXmpi_ZnS ../examples/example00
cp ./mlcoalsimX ../examples/example01
echo cp ./mlcoalsimX ../examples/example01
cp ./mlcoalsimX_ZnS ../examples/example01
echo cp ./mlcoalsimXmpi ../examples/example01
cp ./mlcoalsimXmpi ../examples/example01
echo cp ./mlcoalsimXmpi_ZnS ../examples/example01
cp ./mlcoalsimXmpi_ZnS ../examples/example01
echo cp ./mlcoalsimX ../examples/example10
cp ./mlcoalsimX ../examples/example10
echo cp ./mlcoalsimX_ZnS ../examples/example10
cp ./mlcoalsimX_ZnS ../examples/example10
echo cp ./mlcoalsimXmpi ../examples/example10
cp ./mlcoalsimXmpi ../examples/example10
echo cp ./mlcoalsimXmpi_ZnS ../examples/example10
cp ./mlcoalsimXmpi_ZnS ../examples/example10
echo cp ./mlcoalsimX ../examples/example02
cp ./mlcoalsimX ../examples/example02
echo cp ./mlcoalsimX_ZnS ../examples/example02
cp ./mlcoalsimX_ZnS ../examples/example02
echo cp ./mlcoalsimXmpi ../examples/example02
cp ./mlcoalsimXmpi ../examples/example02
echo cp ./mlcoalsimXmpi_ZnS ../examples/example02
cp ./mlcoalsimXmpi_ZnS ../examples/example02
echo cd ../examples/example00
cd ../examples/example00
echo "sh ./run_mlcoalsim00.sh &"
sh ./run_mlcoalsim00.sh &
echo cd ../example01
cd ../example01
echo "sh ./run_mlcoalsim01.sh &"
sh ./run_mlcoalsim01.sh &
echo cd ../example02
cd ../example02
echo "sh ./run_mlcoalsim02.sh &"
sh ./run_mlcoalsim02.sh &
echo cd ../example10
cd ../example10
echo "sh ./run_mlcoalsim10.sh &"
sh ./run_mlcoalsim10.sh &
echo cd ../../
cd ../../
