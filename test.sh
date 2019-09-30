#!/bin/bash -
#title          :test.sh
#description    :This script test output of ml_coalsim is valid as per available baseline data

set -e

# output folder, case, baseline folder
function validate {
    local output=$1
    local case=$2
    local baseline=$3
    local loci=$4

    pushd $output >/dev/null

    if [[ $loci -eq 1 ]]; then
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
    else
        # N.B.: we must avoid comparing the _summary.out file. It has a timestamp, hence we'll never get a match.
        for file in ${case}A_[l,p]*.out
        do
           cmp --silent $file $baseline/$file
           local status=$?
           if [[ $status -ne 0 ]]; then
            echo >&2 "\t\t Percentiles do not match! ($file ; error code: $status)"
            break
           fi
        done
    fi

    popd >/dev/null

    return $status
}

# $1 the case id (e.g., case01)
# $2 the case file (e.g., 10loci_with_priors)
# $3 the number of loci (e.g., 1)
function run {
    local id=$1
    local case=$2
    local loci=$3
    local output=/tmp/$id

    echo >&2 "-- Running the validation for case '$case' --"

    echo >&2 "\tCreating the folder to put the generated data: $output"
    mkdir -p $output >/dev/null

    
    echo >&2 "\tGenerating the data for case $case"
    local project=$(pwd)
    pushd baselines/$id/ >/dev/null
    $project/bin/mlcoalsimX $case.txt $output/${case}A.txt >/dev/null 2>&1
    popd >/dev/null
    
    echo >&2 "\tValidating the generated data for case $case"

    local here=$(pwd)
    validate $output ${case} $here/baselines/$id $loci
    local status=$?

    echo >&2 "\tTear down"
    
    rm -rf $output

    return $status
}

run case01 1locus_1pop_mhit0_n100_S200 1 || exit 1
run case02 1locus_1pop_mhit0_rec100_n100_S200 1 || exit 1
run case03 1locus_1pop_mhit0_recINF_Lmod 1 || exit 1
run case04 1locus_1pop_mhit0_obs 1 || exit 1
run case05 10loci_with_priors 10 || exit 1

echo "\n***All validations OK!!***"
exit 0
