#!/bin/bash

# Downsamples combined per-barcode FASTQ files from a source directory to a fixed target size
# using bbmap's reformat.sh. 
# Input directory: $BASE_DIR, expected to contain files named like
#   barcodeXX_final_combined.fastq.gz
# Output directory: $SUBSAMPLING_BASE_DIR (/results/downsampling_result); each sample is written as
#   barcodeXX_DS_5.2GB.fastq.gz (target ~5.2 billion bases via samplebasestarget).

# Base directory where the samples are located
BASE_DIR="/media/crowfoot2/DATOS/wastewaters_ch/results/combined_barcodes" #from the concatenation script
SUBSAMPLING_BASE_DIR="/media/crowfoot2/DATOS/wastewaters_ch/results/downsampling_result"

# Create output directory if it doesn't exist
mkdir -p $SUBSAMPLING_BASE_DIR

# Change to bbmap directory
cd /home/crowfoot2/bbmap
chmod +x reformat.sh

for FASTQ_GZ in ${BASE_DIR}/barcode*_final_combined.fastq.gz; do
    # Extract sample name (e.g., barcode01) from the file name
    SAMPLE_NAME=$(basename "$FASTQ_GZ" "_final_combined.fastq.gz")
    
    # Output file path
    OUTPUT_FASTQ_GZ="${SUBSAMPLING_BASE_DIR}/${SAMPLE_NAME}_DS_4.7GB.fastq.gz"
    
    # Check if input FASTQ file exists
    if [[ -f $FASTQ_GZ ]]; then
        # Run reformat.sh to downsample to approximately 5.2 GB
        ./reformat.sh in="$FASTQ_GZ" out="$OUTPUT_FASTQ_GZ" samplebasestarget=5200000000 overwrite=true
    else
        echo "WARNING: ${FASTQ_GZ} does not exist. Skipping ${SAMPLE_NAME}..."
    fi
done

echo "Subsampling processing completed for all samples."
