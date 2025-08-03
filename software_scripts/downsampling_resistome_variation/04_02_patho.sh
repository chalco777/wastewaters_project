#!/bin/bash

# Runs RGI kmer_query on per-barcode RGI main JSON outputs to perform k-mer–
# based pathogen-of-origin prediction for AMR genes.
# Input: BASE_DIR/DS_*/barcode*.json (RGI main output JSON files in results/DS_rgi_result).
# Output: For each barcode, a subdirectory under MAIN_OUTPUT_DIR/DS_<...>/<barcode>/ 
#   containing the kmer_query results (pathogen-of-origin predictions) using the local
#   CARD k-mer database. results/rgi_PATHO_result
# Uses k-mer size 61 (default), requires at least 10 k-mer hits to call a signal,
# and runs in local mode with RGI integration.

# Activate Conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate rgi

# Create the main output directory
MAIN_OUTPUT_DIR="/media/crowfoot2/DATOS/wastewaters_ch/results/rgi_PATHO_result"
mkdir -p $MAIN_OUTPUT_DIR

# Base directory where the DS subdirectories are located
BASE_DIR="/media/crowfoot2/DATOS/wastewaters_ch/results/DS_rgi_result"

 
# Iterate over each DS subdirectory
for DS_DIR in ${BASE_DIR}/DS_*; do
    # Extract DS name (e.g., DS_9) from the directory path
    DS_NAME=$(basename $DS_DIR)
    echo $DS_DIR
    # Create the output directory for the current DS
    DS_OUTPUT_DIR="${MAIN_OUTPUT_DIR}/${DS_NAME}"
    mkdir -p $DS_OUTPUT_DIR
    echo 'here'
    # Iterate over each barcode fastq file in the DS subdirectory
    for JSON in ${DS_DIR}/barcode*.json; do
        # Extract barcode name (e.g., barcode01_9Gb) from the file path
        BARCODE_NAME=$(basename $JSON .json)
        # Create the output directory for the current barcode
        BARCODE_OUTPUT_DIR="${DS_OUTPUT_DIR}/${BARCODE_NAME}"

        # Run RGI on assembled contigs
        rgi kmer_query --rgi --kmer_size 61 --threads 20 --minimum 10 --input $JSON --output $BARCODE_OUTPUT_DIR --local
    done
done

echo "Flye assembly completed for all DS and barcode samples."
