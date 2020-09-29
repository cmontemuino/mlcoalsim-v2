#!/bin/bash -
#title          :validate.sh
#description    :This script test generated files with ml_coalsim-v2 is valid as per available baseline data

set -e

# To make file expansion in for-loop work properly
# https://www.cyberciti.biz/faq/bash-loop-over-file/
# http://www.faqs.org/docs/bashman/bashref_34.html
shopt -s nullglob

# $1 the output folder (e.g., /tmp/example00-123837)
# $2 the case id (e.g., Example1locus_1pop_mhit0)
# $3 the validation folder (e.g., ./validation)
# $4 the case id (e.g., example00)
function validate() {
    local output=$1
    local case=$2
    local validation=$3
    local id=$4

    local baseline="$validation/$id"

    pushd "$output" >/dev/null

    if [[ $id == "example00" ]]; then
        printf >&2 '\t\t3.1 Comparing Percentiles: '
        cmp "${case}_PPercentiles.out" "$baseline/${case}_PPercentiles.out"
        local status=$?
        if [[ $status -eq 0 ]]; then
            printf >&2 '\xE2\x9C\x94\n'
            printf >&2 '\t\t3.2 Comparing Checksums: '
            # Remove first lines containing info about the software version, date, etc.
            gsed -i '1,6d' "${case}.out"
            sha256sum --check --status --strict "$baseline/${case}_SHA256SUMS"
            status=$?
            if [[ $status -ne 0 ]]; then
                printf >&2 '  \xE2\x9C\x98 (error code: %s)\n' "$status"
            else
                printf >&2 '  \xE2\x9C\x94\n'
            fi
        else
            printf >&2 '\xE2\x9C\x98 (error code: %s)\n' "$status"
        fi
    else
        # Remove first lines containing info about the software version, date, etc.
        [[ -f "$case".out ]] && gsed -i '1,6d' "$case".out
        
        for file in "$case"_summary.out; do [[ -f "$file" ]] && gsed -i '1,6d' "$file"; done

        for file in "$case"_*_rank[0-9][0-9][0-9].out; do gsed -i '1,1d' "$file"; done

        # This handles the special case of MPI with locus where "rank" files are not generated.
        for file in "$case"_locus_[0-9][0-9][0-9][0-9][0-9].out; do gsed -i '1,1d' "$file"; done

        # We make a file-to-file comparison
        printf >&2 '\t\t3.1 Comparing Output: '

        for file in "$case"*.out
        do
            cmp --silent "$file" "$baseline/$file"
            local status=$?
            if [[ $status -ne 0 ]]; then
                printf >&2 '\t\t\xE2\x9C\x98 Output does not match! (%s, error code: %s)\n' "$file" "$status"
                break
            fi
        done
        if [[ $status -eq 0 ]]; then printf >&2 '  \xE2\x9C\x94\n'; fi
    fi

    popd >/dev/null

    return $status
}

# $1 the case id (e.g., example00)
# $2 the case file (e.g., Example1locus_1pop_mhit0)
# $3 the variant (e.g., standard, zns, mpi)
# $4 the number of processors to run mpi (i.e., `np`)
function run() {
    local id=$1
    local case=$2
    local variant=$3
    local np=$4

    printf >&2 '** Running the validation for |%s|, variant |%s| **\n' "$case" "$variant"

    # Maintain a reference to the current location that we later use for running the simulations
    # and validate the results.
    local project
    project=$(pwd)

    # Create a temporary folder to put the outputs of running this simulation
    local output="/tmp/$id-$RANDOM"

    # Letter used as suffix in the output files. Usually "A", but "B" in case of ZNS variant.
    # This distinction only applies to example00
    local letter
    if [[ $id == "example00" ]]; then letter="A"; fi

    local status=0
    case $variant in
      zns)
        checkBinary "$project/build/mlcoalsimX_ZnS" && status=0 || status=1
        COMMAND="$project/build/mlcoalsimX_ZnS"
        if [[ $id == "example00" ]]; then letter="B"; fi
        ;;
      mpi)
        checkBinary "$project/build/mlcoalsimXmpi" && status=0 || status=1
        COMMAND="mpirun -np $np $project/build/mlcoalsimXmpi"
        ;;
      *)
        checkBinary "$project/build/mlcoalsimX" && status=0 || status=1
        COMMAND="$project/build/mlcoalsimX"
        ;;
    esac

    if [[ $status -eq 0 ]]; then

      printf >&2 '\t1. Creating the folder to put the generated data: %s\n' "$output"
      mkdir -p "$output" >/dev/null

      printf >&2 '\t2. Generating the data\n'
      pushd "examples/$id/" >/dev/null

      $COMMAND "$case.txt" "$output/${case}${letter}.out" >/dev/null 2>&1

      popd >/dev/null

      printf >&2 '\t3. Validating the generated data for |%s|\n' "$case"

      validate "$output" "${case}${letter}" "$project/validation" "$id"
      status=$?

      printf >&2 "\t4. Tear down\n"

      rm -rf "$output"
    else
      printf >&2 '\t Cannot run |%s|. You need to generate binaries!!\n' "$COMMAND"
    fi

    return $status
}


function checkBinary() {
  local _binary="$1"

  # Checks if the binary is available.
  command -v "$_binary"  >/dev/null
  local status=$?
  return $status
}

OK=0
KO=0
function example00() {
    run example00 Example1locus_1pop_mhit0 standard && OK=$((OK + 1)) || KO=$((KO + 1))
    run example00 Example1locus_1pop_mhit0_n100_S200 standard && OK=$((OK + 1)) || KO=$((KO + 1))
    run example00 Example1locus_1pop_mhit0_rec100_n100_S200 standard && OK=$((OK + 1)) || KO=$((KO + 1))
    run example00 Example1locus_1pop_mhit0_obs zns && OK=$((OK + 1)) || KO=$((KO + 1))
    run example00 Example1locus_1pop_mhit0_rec100_n100_S200 standard && OK=$((OK + 1)) || KO=$((KO + 1))
    run example00 Example1locus_1pop_mhit0_rec100_t10_n100 standard && OK=$((OK + 1)) || KO=$((KO + 1))
    run example00 Example1locus_1pop_mhit0_recINF_Lmod standard && OK=$((OK + 1)) || KO=$((KO + 1))
    run example00 Example1locus_1pop_mhit0_recINF_S200_Lfix standard && OK=$((OK + 1)) || KO=$((KO + 1))
}

function example01() {
    run example01 Example1locus standard && OK=$((OK + 1)) || KO=$((KO + 1))
}

function example02() {
    run example02 mlcoal_input1 zns && OK=$((OK + 1)) || KO=$((KO + 1))
    run example02 mlcoal_input2 zns && OK=$((OK + 1)) || KO=$((KO + 1))
}

function example10() {
    run example10 Example10loci mpi 4 && OK=$((OK + 1)) || KO=$((KO + 1))
}

example00
example01
example02
example10

printf >&2 '\n***Passed: %d***\n' "$OK"
printf >&2 '***Not Passed: %d***\n' "$KO"
exit 0
