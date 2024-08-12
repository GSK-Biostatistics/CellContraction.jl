#!/bin/bash

# Default values
julia_exec=""
project_path=""
julia_script=""
nworkers=""
template_file=""

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --template=*)
            template_file="${key#*=}"
            shift
            ;;
        --julia=*)
            julia_exec="${key#*=}"
            shift
            ;;
        --project=*)
            project_path="${key#*=}"
            shift
            ;;
        --script=*)
            julia_script="${key#*=}"
            shift
            ;;
        --nworkers=*)
            nworkers="${key#*=}"
            shift
            ;;
        --idx=*)
            cmpd_mech_idx="${key#*=}"
            shift
            ;;
        *)
            echo "Unknown option: $key"
            exit 1
            ;;
    esac
done

# Check if required arguments are provided
if [[ -z "$template_file" || -z "$julia_exec" || -z "$project_path" || -z "$julia_script" || -z "$nworkers" ]]; then
    echo "Usage: ./generate_slurm_from_template.sh --template=<template_file_path> --julia=<julia_binary_path> --project=<project_path> --script=<julia_script_path> --nworkers=<num_workers> [--idx=<compound_or_mechanism_index_range>]"
    exit 1
fi

# Read the SLURM template
template=$(<$template_file)

# Replace placeholders with arguments
template=${template//"\${JULIA}"/$julia_exec}
template=${template//"\${PROJECT}"/$project_path}
template=${template//"\${SCRIPT}"/$julia_script}
template=${template//"\${NWORKERS}"/$nworkers}
template=${template//"\${IDX}"/$cmpd_mech_idx}

# Save the generated SLURM script to a file
filename=$(basename "$julia_script" .jl)
output_file="${filename}.slurm"
echo "$template" > "$output_file"

echo "Generated SLURM script saved to $output_file"
