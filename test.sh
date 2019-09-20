#!/bin/bash -
#title          :test.sh
#description    :This script test output of ml_coalsim is valid as per available baseline data

set -e

# output folder, case, baseline folder
function validate {
    local output=$1
    local case=$2
    local baseline=$3

    pushd $output >/dev/null

    cmp --silent ${case}A_PPercentiles.out $baseline/${case}A_PPercentiles.out
    local status=$?
    if [[ $status -eq 0 ]]; then
        sha256sum --check --status --strict $baseline/${case}A_linkedlocus_sha256.txt
        status=$?
        if [[ $status -ne 0 ]]; then
            echo >&2 "\t\t Checksum do not match! (error code: $status)"
        fi
    else
        echo >&2 "\t\t Percentiles do not match! (error code: $status)"
    fi

    popd >/dev/null

    return $status
}

# the case number
function run {
    local id=$1
    local case=$2
    local output=/tmp/$id

    echo >&2 "-- Running the validation for case '$case' --"

    mkdir -p $output

    echo >&2 "\tCreating the folder to put the generated data: $output"
    bin/mlcoalsimX baselines/$id/$case.txt $output/${case}A.txt >/dev/null 2>&1

    echo >&2 "\tValidating the generated data for case $case"
    local here=$(pwd)
    validate $output ${case} $here/baselines/$id
    local status=$?

    echo >&2 "\tTear down"
    
    rm -rf $output

    return $status
}

run case01 1locus_1pop_mhit0_n100_S200 || exit 1
run case02 1locus_1pop_mhit0_rec100_n100_S200 || exit 1
run case03 1locus_1pop_mhit0_recINF_Lmod || exit 1
run case04 1locus_1pop_mhit0_obs || exit 1

echo "\n***All validations OK!!***"
exit 0
