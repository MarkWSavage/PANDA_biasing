#!/bin/bash
set -uo pipefail
cd /home/markr/geant4/PANDA_biasing

# Wait for the already-running neutron biased job (started before this
# script existed) to finish before touching Results/Current/ at all --
# only one PANDA process may write to that hardcoded path at a time.
NEUTRON_BIASED_PID=8880
echo "=== WAITING for in-flight neutron biased run (PID $NEUTRON_BIASED_PID) ==="
while kill -0 "$NEUTRON_BIASED_PID" 2>/dev/null; do
    sleep 2
done
echo "=== neutron biased run process exited ==="

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

# Archive the neutron biased run that just finished (its raw stdout/log
# already exists from the earlier launch).
mkdir -p Results/Neutron_1041_validated
mv Results/Current/events.csv Results/Neutron_1041_validated/events.csv
grep -A6 "GLOBAL FINAL COUNTS" run_output_neutron_biased.log > Results/Neutron_1041_validated/summary.txt
echo "=== ARCHIVED: Results/Neutron_1041_validated/ ==="

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
