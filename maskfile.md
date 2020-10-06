# Tasks

Please notice this is a [mask markdown file][mask]. You need nothing special for reading it
and using the commands provided here, but if you want a bit more of automation,
then you will need to [install mask][mask-install].

If you choose to do so, then running `mask --help` will give you an overview of the tasks you can run.

## build

> Builds the application with `CMAKE` and genearate the executables. An optional `--target` flag can be passed to specify a single target.

It is possible to tweak the flags used by `CMAKE` by specifying the `btype` flag:

* `Debug`: Adds the `-g` flag
* `Release`: Adds the `-O3 -DNDEBUG` flags to the compiler
* `MinSizeRel`: Adds `-Os -DNDEBUG` flags
* `RelWithDebInfo`: Adds `-O2 -g -DNDEBUG` flags (**default value**)

Please note that all executables depend on `zf_logs` and `benchmark`, therefore, even if you issue `mask build -t mlcoalsimXmpi`,
then you will notice the following:

```sh
mask build -t mlcoalsimXmpi
-- Configuring done
-- Generating done
-- Build files have been written to: <some_path_to_mlcoalsim>/build
Building target <mlcoalsimXmpi> ...
[ 15%] Built target zf_log
[ 30%] Built target benchmark
[100%] Built target mlcoalsimXmpi
```

**OPTIONS**
* target
  * flags: -t --target
  * type: string
  * desc: Which target to build (valid targets: `mlcoalsimX`, `mlcoalsimX_ZnS`, `mlcoalsimXmpi`, `mlcoalsimXmpi_ZnS`, `mlcoalsimB`, `mlcoalsimBmpi`)
* btype
  * flags: -b --btype
  * type: string
  * desc: Build type used by CMAKE (valid types: `Debug`, `Release`, `MinSizeRel`, `RelWithDebInfo`). Default is `RelWithDebInfo`. 

~~~sh

_target=${target:-all}
_build_type=${btype:-RelWithDebInfo}

mkdir -p build

cmake -B build -DCMAKE_BUILD_TYPE="$_build_type"

if [[ "$_target" == "all" ]]; then
    echo "Building all targets..."
    make -j4 -C build
else
    printf 'Building target <%s> ...\n' "$_target"
    make -j4 -C build "$_target"
fi
~~~

## validate

> Validates the application works as expected from a functional point of view.
> There are several subcommands that you might discover by issuing `mask validate --help`.  

Source code changes are common as part of the lifecycle of whatever project. As is in the case of many
bioinformatics tools, unit and integration tests are very limited. This project is no different.

Source code form this project considerably diverged from the [original mlcoalsimv2][mlcoalsimv2-original]. Therefore, we need to verify the
 results are still valid.

For percentiles files we can directly compare them with the [cmp Linux tool][cmp], but for output files we use the [sha256sum Linux tool][sha256sum].

We have ran some examples provided in the [examples](examples) folder, and generated a set of checksums and
output files that you can find in the [validation](validation) folder.

We provide validation scripts to make sure that the application behaves as expected.

**Generating the Checksum**

We needed to massage the output files a bit before creating the SHA256 checksum. More specifically, we need to remove the first 6 lines that contain
data related to the version of the program, date when the file was created, and the input file. Then we can create the checksum:

_**Important Note for Mac OSX users**_: I have GNU sed installed and linked to `gsed`. If don't have it installed, then you will need to provide the
extension for backup files when using the "in-place replacement" feature. See more details [here][sed-macosx-issue].

```shell script
build/mlcoalsimX examples/example00/Example1locus_1pop_mhit0.txt build/Example1locus_1pop_mhit0A.txt
pushd build
gsed -i '1,6d' Example1locus_1pop_mhit0A.txt
sha256sum Example1locus_1pop_mhit0A.txt | tee ../validation/example00/Example1locus_1pop_mhit0A_SHA256SUMS
mv Example1locus_1pop_mhit0A_PPercentiles.out ../validation/example00/.
popd
```

**Generate validations for example01**

In the case of the `example01` data, that is using "prior" files, we need to slightly change what we have been done so far.

```shell script
pushd examples/example01
../../build/mlcoalsimX ./Example1locus.txt Example1locus.out
gsed -i '1,6d' Example1locus_summary.out
gsed -i '1,1d' Example1locus_locus_00000_rank000.out
mv Example1locus_*.out ../../validation/example01/.
popd
```

**Generate validations for example10**

In the case of `example10` data, we have "prior" files and MPI. The generation of validation files is different here too:

```shell script
pushd examples/example10
mpirun -np 4 ../../build/mlcoalsimXmpi Example10loci.txt Example10loci.out
gsed -i '1,6d' Example10loci_summary.out
for file in $(ls Example10loci_locus_*.out); do gsed -i '1,1d' $file; done
mv Example10loci_*.out ../../validation/example10/.
popd
```

**Validate the Output is as Expected**

Now we need to generate the outputs with the new executable and perform the verifications:

```shell script
build/mlcoalsimX examples/example00/Example1locus_1pop_mhit0_n100_S200.txt build/Example1locus_1pop_mhit0_n100_S200A.txt
pushd build
gsed -i '1,6d' Example1locus_1pop_mhit0_n100_S200A.txt
grep Example1locus_1pop_mhit0_n100_S200A.txt ../validation/example00/Example1locus_1pop_mhit0_n100_S200A_SHA256SUMS | tee /dev/fd/2 | sha256sum --check --strict  -
# 3176b276d45244dbe4f98760b1ca3fdd939e0580bd82198848c495d86471b603  Example1locus_1pop_mhit0A.txt
# Example1locus_1pop_mhit0A.txt: OK

cmp Example1locus_1pop_mhit0_n100_S200A_PPercentiles.out ../validation/example00/Example1locus_1pop_mhit0_n100_S200A_PPercentiles.out
# Output should be empty!!!
popd
```

**Validate simulations from example01 and example02**

We are not using checksum files for `example01` and `example02`, hence comparison with the `cmp` command is enough:

```shell script
build/mlcoalsimX_ZnS build/mlcoalsimX_ZnS examples/example02/mlcoal_input1 build/mlcoal_output1.txt
build/mlcoalsimX_ZnS build/mlcoalsimX_ZnS examples/example02/mlcoal_input2 build/mlcoal_output2.txt
pushd build
gsed -i '1,6d' mlcoal_output1.txt
gsed -i '1,6d' mlcoal_output2.txt
gsed -i '1,1d' mlcoal_output2_linkedlocus__rank000.out

cmp mlcoal_output1.txt ../validation/example02/mlcoal_output1.txt
# Output should be empty!!!
cmp mlcoal_output1_PPercentiles.out ../validation/example02/mlcoal_output1_PPercentiles.out
# Output should be empty!!!

cmp mlcoal_output2.txt ../validation/example02/mlcoal_output2.txt
# Output should be empty!!!
cmp mlcoal_output2_linkedlocus__rank000.out ../validation/example02/mlcoal_output2_linkedlocus__rank000.out
# Output should be empty!!!
popd
```

**Validate simulations from example10**

```shell script
pushd examples/example10
mpirun -np 4 ../../build/mlcoalsimXmpi Example10loci.txt Example10loci.out
gsed -i '1,6d' Example10loci_summary.out
for file in $(ls Example10loci_locus_*.out); do gsed -i '1,1d' $file; done

for file in $(ls Example10loci_*.out); do echo "Validating |$file|" && cmp $file ../../validation/example10/$file; done
cmp Example10loci_summary.out ../../validation/example10/Example10loci_summary.out
# Output should be empty!!!
rm -f Example10loci_*.out
popd
```

### all

> Run all the validations. This task might take long, but it might be the best way to validate
> the expected functionality.

~~~sh
validation/validate.sh
~~~

[mask]: https://github.com/jakedeichert/mask
[mask-install]: https://github.com/jakedeichert/mask/#installation

