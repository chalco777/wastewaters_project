#!/bin/bash

# Runs RGI main on Flye assemblies to predict antibiotic resistance genes.
# Input: per-barcode assembly.fasta files under BASE_DIR (assemblies produced by Flye).
# Output: RGI results (JSON/text) per barcode under MAIN_OUTPUT_DIR/<barcode>/ using the
#        local CARD database. Uses 28 threads and the default detection settings
#        (Perfect and Strict hits, with nudging enabled, no loose hits unless promoted).
# Requires: an accessible localDB/ directory (RGI --local), RGI installed in the activated
#           conda environment, and that the script is executed from the correct parent of localDB/.


# Activate Conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate rgi


# Create the main output directory
MAIN_OUTPUT_DIR="/home/crowfoot2/DATOS/wastewaters_ch/results/DS_rgi_result"
mkdir -p $MAIN_OUTPUT_DIR

# Base directory where the DS subdirectories are located
BASE_DIR="/media/crowfoot2/DATOS/wastewaters_ch/results/DS_flye"

 
# Iterate over each DS subdirectory
for DS_DIR in $BASE_DIR/barcode*; do
    # Extract DS name (e.g., DS_9) from the directory path
    DS_NAME=$(basename $DS_DIR)
    echo "Procesando $DS_DIR"
    # Create the output directory for the current DS
    DS_OUTPUT_DIR="$MAIN_OUTPUT_DIR/$DS_NAME"
    mkdir -p $DS_OUTPUT_DIR
    echo 'here'
    ASSEMBLY_FASTA="$DS_DIR/assembly.fasta"
    rgi main --input_sequence $ASSEMBLY_FASTA --output_file $DS_OUTPUT_DIR --local -n 28 
done


echo "Flye assembly completed for all DS and barcode samples."
