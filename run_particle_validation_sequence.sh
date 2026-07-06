#!/bin/bash
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

run_and_archive() {
    local mac="$1"
    local archive_dir="$2"
    local log="$3"

    echo "=== STARTING: $archive_dir ($mac) ==="
    ./PANDA "$mac" > "$log" 2>&1
    local rc=$?
    echo "=== FINISHED: $archive_dir (exit $rc) ==="

    if [ $rc -ne 0 ]; then
        echo "=== ERROR: $archive_dir exited non-zero, stopping sequence ==="
        exit 1
    fi

    mkdir -p "Results/$archive_dir"
    mv Results/Current/events.csv "Results/$archive_dir/events.csv"
    grep -A6 "GLOBAL FINAL COUNTS" "$log" > "Results/$archive_dir/summary.txt"
    echo "=== ARCHIVED: Results/$archive_dir/ ==="
}

run_and_archive Macros/run_neutron_biased.mac   Neutron_1041_validated      run_output_neutron_biased.log
run_and_archive Macros/run_neutron_unbiased.mac Neutron_unbiased_validated run_output_neutron_unbiased.log

run_and_archive Macros/run_alpha_biased.mac      Alpha_1041_validated      run_output_alpha_biased.log
run_and_archive Macros/run_alpha_unbiased.mac    Alpha_unbiased_validated  run_output_alpha_unbiased.log

run_and_archive Macros/run_deuteron_biased.mac   Deuteron_1041_validated     run_output_deuteron_biased.log
run_and_archive Macros/run_deuteron_unbiased.mac Deuteron_unbiased_validated run_output_deuteron_unbiased.log

run_and_archive Macros/run_triton_biased.mac     Triton_1041_validated     run_output_triton_biased.log
run_and_archive Macros/run_triton_unbiased.mac   Triton_unbiased_validated run_output_triton_unbiased.log

run_and_archive Macros/run_He3_biased.mac        He3_1041_validated     run_output_He3_biased.log
run_and_archive Macros/run_He3_unbiased.mac      He3_unbiased_validated run_output_He3_unbiased.log

echo "=== ALL_RUNS_COMPLETE ==="
