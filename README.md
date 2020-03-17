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

## How to Build

First we need to build the project. In the following script we are creating a `build` folder in
the root of this project (which by the way is already git-ignored). Please feel free to
use whatever other folder if you like to:

```shell script
mkdir build
cmake -B build -DCMAKE_BUILD_TYPE=<build_type>
make -C build
```

where `<build_type>` must be some of:

* `Debug`: Adds the `-g` flag
* `Release`: Adds the `-O3 -DNDEBUG` flags to the compiler
* `MinSizeRel`: Adds `-Os -DNDEBUG` flags
* `RelWithDebInfo`: Adds `-O2 -g -DNDEBUG` flags  **_<--- Default option if `DCMAKE_BUILD_TYPE` is not provided_**

Four executables are going to be generated:
* `mlcoalsimX`: application running without MPI
* `mlcoalsimX_ZnS`: application running without MPI but including ZnS statistics
* `mlcoalsimXmpi`: application running with MPI
* `mlcoalsimXmpi_ZnS`: application running with MPI and including ZnS statistics

If you want to generate one executable only, then just reference it when calling `make`:

```shell script
make -C build mlcoalsimX_ZnS
```

## How to Run the Examples

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

Other examples you might want to run:

* `build/mlcoalsimX examples/example00/Example1locus_1pop_mhit0_rec100.txt build/Example1locus_1pop_mhit0_rec100.out`
* `mpirun -np 4 build/mlcoalsimXmpi examples/example00/Example1locus_1pop_mhit0_rec100.txt build/Example1locus_1pop_mhit0_rec100.out`
* `build/mlcoalsimX_ZnS examples/example00/Example1locus_1pop_mhit0_rec100_S20_n100.txt build/Example1locus_1pop_mhit0_rec100_S20_n100.out`
* `mpirun -np 4 build/mlcoalsimXmpi_ZnS examples/example00/Example1locus_1pop_mhit0_rec100_S20_n100.txt build/Example1locus_1pop_mhit0_rec100_S20_n100.out`

## How to Verify the Results

Source code form this project considerably diverged from the [original mlcoalsimv2][mlcoalsimv2-original]. Therefore, we need to
verify the results are still valid.

For percentiles files we can directly compare them with the [cmp Linux tool][cmp], but for output files we use the [sha256sum Linux tool][sha256sum].
Output files can grow very large (in the order of hundreds of megabytes) depending on the .

### Generating the Checksum

We need to massage the output files a bit before creating the SHA256 checksum. More specifically, we need to remove the first 6 lines that contain
data related to the version of the program, date when the file was created, and the input file. Then we can create the checksum:

**Important Note for Mac OSX users**: I have GNU sed installed and linked to `gsed`. If don't have it installed, then you will need to provide the
extension for backup files when using the "in-place replacement" feature. See more details [here][sed-macosx-issue].

```shell script
build/mlcoalsimX examples/example00/Example1locus_1pop_mhit0.txt build/Example1locus_1pop_mhit0A.txt
pushd build
gsed -i '1,6d' Example1locus_1pop_mhit0A.txt
sha256sum Example1locus_1pop_mhit0A.txt | tee ../validation/Example1locus_1pop_mhit0A_SHA256SUMS
mv Example1locus_1pop_mhit0A_PPercentiles.out ../validation/.
popd 
```

### Validate the Output is as Expected

Now we need to generate the outputs with the new executable and perform the verifications:

```shell script
build/mlcoalsimX examples/example00/Example1locus_1pop_mhit0.txt build/Example1locus_1pop_mhit0A.txt
pushd build
gsed -i '1,6d' Example1locus_1pop_mhit0A.txt
grep Example1locus_1pop_mhit0A.txt ../validation/Example1locus_1pop_mhit0A_SHA256SUMS | tee /dev/fd/2 | sha256sum --check --strict  -
# 3176b276d45244dbe4f98760b1ca3fdd939e0580bd82198848c495d86471b603  Example1locus_1pop_mhit0A.txt
# Example1locus_1pop_mhit0A.txt: OK

cmp Example1locus_1pop_mhit0A_PPercentiles.out ../validation/Example1locus_1pop_mhit0A_PPercentiles.out
# Output should be empty!!!
popd 
```



[clion]: https://www.jetbrains.com/clion/
[cmake]: https://cmake.org
[cmp]: https://linux.die.net/man/1/cmp
[homebrew]: https://brew.sh
[homebrew-custom-openmpi]: https://github.com/cmontemuino/homebrew-custom
[mlcoalsimv2-original]: https://github.com/CRAGENOMICA/mlcoalsim-v2
[openmpi]: https://www.open-mpi.org
[sed-macosx-issue]: https://stackoverflow.com/questions/4247068/sed-command-with-i-option-failing-on-mac-but-works-on-linux/4247319#4247319
[sha256sum]: https://linux.die.net/man/1/sha256sum/