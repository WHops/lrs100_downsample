#!/bin/bash

# Usage function
usage() { echo "Usage: $0 -t <threads>"; exit 1; }

# Parse threads option
while getopts "t:" opt; do
    case $opt in
        t) threads=$OPTARG ;;
        *) usage ;;
    esac
done

# Ensure threads parameter is provided
[ -z "$threads" ] && usage

# Hard-coded parameters
data_basedir="/ifs/data/research/projects/jordi/wp130_downsample_Phase1_Revio/results/permutation"
X_values=("10X" "15X" "20X")
permutations=({0..9})


# Generate commands and run them in parallel
commands=()
for X in "${X_values[@]}"; do
    for perm in "${permutations[@]}"; do
        output="variants_report_${X}_${perm}.tsv"
        commands+=("python scripts/find_lrs100_variants.py data/fullvar_withvar.tsv ${data_basedir}/${X} ${output} ${perm}")
    done
done

# Run the commands with the specified number of parallel threads
printf "%s\n" "${commands[@]}" | xargs -n 1 -P "$threads" -I {} bash -c '{}'
