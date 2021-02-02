#!/bin/bash

function configure() {
    cat <<EOF
    {
    "input": "wkt/prague.wkt",
    "pers_beta":5.0,
    "dp_eps":5.0,
    "sds_beta":5.0,
    "sds_eps":2.5,
    "sds_rounds":5,
    "mrs_beta":5.0,
    "mrs_eps":2.5,
    "mrs_rounds":5.0,
    "mrs_rounds":5
}

EOF
}

export -f configure

./main <(configure) profile.json --profile
