# mlcoalsim-v2

## Quick Start

This project requires both [OpenMPI][openmpi] and [CMake][cmake].

If you are running on Mac OSX and you happen to use [Homebrew][homebrew], then you might want to check
a custom for installing OpenMPI [here][homebrew-custom-openmpi].

## Enabling MPI in the IDE

If you open the project with a contextual IDE, for example [CLion][clion], then you will
notice that all the code related to MPI will not be "clickable". One way to resolve this
is by setting the environment variable `WITH_MPI`. You just need to provide whatever value to it.

*Note: what I do is setting the environment variable inside the IDE, so that I do not have problems
with building the project (Preferences -> Build, Execution, Deployment -> CMake -> Environment: WITH_MPI=1 ZNS_ACTIVE=1).*

## How to Build

Please refer to the [maskfile.md file](maskfile.md) that contains a section for [building the project](maskfile.md#build).

## How to Run the Examples

Several examples are provided in the `examples` folder. You can run most of them in the following way:

```shell script
mpirun -np 4 build/mlcoalsimXmpi_ZnS  examples/example00/Example1locus_1pop_mhit0_rec100.txt build/Example1locus_1pop_mhit0_rec100.out

# Without MPI:
# build/mlcoalsimX  examples/example00/Example1locus_1pop_mhit0_rec100.txt build/Example1locus_1pop_mhit0_rec100.out
```

In the case of the `example01` and `example10`, where a "prior" file is being used, you need to switch the directory first.
For example:

```shell script
pushd examples/example10
mpirun -np 4 ../../build/mlcoalsimXmpi  Example10loci.txt ../../build/Example10loci.out
popd
```

Other examples you might want to run:

* `build/mlcoalsimX examples/example00/Example1locus_1pop_mhit0_rec100.txt build/Example1locus_1pop_mhit0_rec100.out`
* `mpirun -np 4 build/mlcoalsimXmpi examples/example00/Example1locus_1pop_mhit0_rec100.txt build/Example1locus_1pop_mhit0_rec100.out`
* `build/mlcoalsimX_ZnS examples/example00/Example1locus_1pop_mhit0_rec100_S20_n100.txt build/Example1locus_1pop_mhit0_rec100_S20_n100.out`
* `mpirun -np 4 build/mlcoalsimXmpi_ZnS examples/example00/Example1locus_1pop_mhit0_rec100_S20_n100.txt build/Example1locus_1pop_mhit0_rec100_S20_n100.out`

## How to Verify the Results

Please refer to the [validation section in maskfile.md file](maskfile.md#validate).


[clion]: https://www.jetbrains.com/clion/
[cmake]: https://cmake.org
[cmp]: https://linux.die.net/man/1/cmp
[homebrew]: https://brew.sh
[homebrew-custom-openmpi]: https://github.com/cmontemuino/homebrew-custom
[mlcoalsimv2-original]: https://github.com/CRAGENOMICA/mlcoalsim-v2
[openmpi]: https://www.open-mpi.org
[sed-macosx-issue]: https://stackoverflow.com/questions/4247068/sed-command-with-i-option-failing-on-mac-but-works-on-linux/4247319#4247319
[sha256sum]: https://linux.die.net/man/1/sha256sum/
