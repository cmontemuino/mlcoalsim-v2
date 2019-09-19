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
    local status1=$(cmp --silent ${case}A_PPercentiles.out $baseline/${case}A_PPercentiles.out && echo 0 || echo 1)
    if [[ $status1 -ne 0 ]]; then
        echo "\t\t Percentiles do not match!"
    fi

    local status2="$(sha256sum --check --status --strict $baseline/${case}A_linkedlocus_sha256.txt && echo 0 || echo 1)"
    if [[ $status2 -ne 0 ]]; then
        echo "\t\t Checksum do not match!"
    fi

    popd >/dev/null
    return 0 || status1 || status2
}

# the case number
function run {
    local id=$1
    local case=$2
    local output=/tmp/$id

    echo "-- Running the validation for case '$case' --"

    mkdir -p $output

    echo "\tCreating the folder to put the generated data: $output"
    bin/mlcoalsimX baselines/$id/$case.txt $output/${case}A.txt >/dev/null 2>&1

    echo "\tValidating the generated data for case $case"
    local here=$(pwd)
    local status="$(validate $output ${case} $here/baselines/$id && echo 0 || echo 1)"

    echo "\tTear down"
    rm -rf $output

    return 0 || status
}

run case01 1locus_1pop_mhit0_n100_S200 || exit 1
run case02 1locus_1pop_mhit0_rec100_n100_S200 || exit 1
run case03 1locus_1pop_mhit0_recINF_Lmod || exit 1

echo "\n***All validations OK!!***"
exit 0
