# mlcoalsim-v2

# Quick Start

This project requires both [OpenMPI][openmpi] and [CMake][cmake].

If you are running on Mac OSX and you happen to use [Homebrew][homebrew], then you might want to check
a custom for installing OpenMPI [here][homebrew-custom-openmpi].

## Enabling MPI in the IDE

If you open the project with a contextual IDE, for example [CLion][clion], then you will
notice that all the code related to MPI will not be "clickable". One way to resolve this
is by setting the environment variable `WITH_MPI`. You just need to provide whatever value to it.

*Note: what I do is setting the environment variable inside the IDE, so that I do not have problems
with building the project.* 

## Building

First we need to build the project. In the following script we are creating a `build` folder in
the root of this project (which by the way is already git-ignored). Please feel free to
use whatever other folder if you like to:

```shell script
mkdir build
cmake -B build
make -C build
```

Four executables are going to be generated:
* `mlcoalsimX`: application running without MPI
* `mlcoalsimX_ZnS`: application running without MPI but including ZnS statistics
* `mlcoalsimXmpi`: application running with MPI
* `mlcoalsimXmpi_ZnS`: application running with MPI and including ZnS statistics

If you want to generate one executable only, then just reference it when calling `make`:

```shell script
make -C build mlcoalsimX_ZnS
```

## Running the Examples

Several examples are provided in the `examples` folder. You can run most of them in the following way:

```shell script
mpirun -np 4 build/mlcoalsimXmpi_ZnS  examples/example00/Example1locus_1pop_mhit0_rec100.txt build/Example1locus_1pop_mhit0_rec100.out

# Without MPI:
# build/mlcoalsimX  examples/example00/Example1locus_1pop_mhit0_rec100.txt build/Example1locus_1pop_mhit0_rec100.out
```

In the case of the `example10`, which uses a "prior" file, then you'll need to switch the directory first.
For example:

```shell script
pushd examples/example10
mpirun -np 4 ../../build/mlcoalsimXmpi  Example10loci.txt ../../build/Example10loci.out
popd
```

[clion]: https://www.jetbrains.com/clion/
[cmake]: https://cmake.org
[homebrew]: https://brew.sh
[homebrew-custom-openmpi]: https://github.com/cmontemuino/homebrew-custom
[openmpi]: https://www.open-mpi.org