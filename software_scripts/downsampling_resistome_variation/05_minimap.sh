#!/bin/bash

# Aligns downsampled nanopore reads back to their corresponding Flye assemblies,
# producing sorted and indexed BAMs per barcode.
# Input:
#   - BASE_DIR: directory with barcode FASTQ files (e.g., downsampled nanopore reads)
#       named like barcodeXXX*.fastq.gz. i.e. /results/downsampling_result
#   - ASSEMBLY_DIR: directory with per-barcode Flye assembly subdirectories
#       (matching the barcode names) each containing assembly.fasta. i.e. /results/DS_flye
# Output:
#   - OUTPUT_DIR/<barcode>/: contains <barcode>_sorted.bam and its index (from samtools). i.e. /results/minimap_result
# The script activates the conda environment (assumed to provide samtools/minimap2),
# runs minimap2 with the map-ont preset, pipes into samtools to convert to BAM,
# sorts and indexes the alignments.


# Define the base directories
BASE_DIR="/media/crowfoot2/DATOS/wastewaters_ch/results/downsampling_result"
ASSEMBLY_DIR="/media/crowfoot2/DATOS/wastewaters_ch/results/DS_flye"
OUTPUT_DIR="/media/crowfoot2/DATOS/MW_P2/wastewaters_ch/results/minimap_result"

# Check if the required directories exist
if [[ ! -d $BASE_DIR || ! -d $ASSEMBLY_DIR ]]; then
    echo "ERROR: One or more required directories do not exist."
    exit 1
fi

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate medaka

# Iterate over each barcode fastq file in the base directory
for FASTQ_FILE in "${BASE_DIR}"/barcode*.fastq.gz; do
    # Check if FASTQ_FILE exists
    if [[ ! -f "$FASTQ_FILE" ]]; then
        echo "No FASTQ files found in ${BASE_DIR}."
        exit 1
    fi

    # Extract barcode name (e.g., barcode18_1405Gb) from the file path
    BARCODE_NAME=$(basename "$FASTQ_FILE" .fastq.gz)
    
    # Path for the corresponding assembly.fasta
    ASSEMBLY_PATH="${ASSEMBLY_DIR}/${BARCODE_NAME}/assembly.fasta"
    
    # Check if the assembly.fasta file exists
    if [[ -f "$FASTQ_FILE" && -f "$ASSEMBLY_PATH" ]]; then
        # Ensure the minimap2 result directory for the sample exists
        mkdir -p "${OUTPUT_DIR}/${BARCODE_NAME}"
        
        # Log the start of processing
        echo "Processing ${BARCODE_NAME}..."

        # Path for the output BAM and sorted BAM files
        BAM_FILE="${OUTPUT_DIR}/${BARCODE_NAME}/${BARCODE_NAME}.bam"
        SORTED_BAM_FILE="${OUTPUT_DIR}/${BARCODE_NAME}/${BARCODE_NAME}_sorted.bam"

        # Run Minimap2 and pipe directly to samtools to output BAM, sort, and index
        minimap2 -ax map-ont -t 30 "$ASSEMBLY_PATH" "$FASTQ_FILE" | \
        samtools view -Sb - | \
        samtools sort -o "$SORTED_BAM_FILE"
        if [[ $? -ne 0 ]]; then
            echo "ERROR: Minimap2 or samtools failed for ${BARCODE_NAME}. Check the input files."
            continue
        fi

        # Index the sorted BAM file
        samtools index "$SORTED_BAM_FILE"
        if [[ $? -ne 0 ]]; then
            echo "ERROR: BAM indexing failed for ${BARCODE_NAME}."
            continue
        fi

        echo "Completed processing ${BARCODE_NAME}."
    else
        echo "WARNING: Either ${FASTQ_FILE} or ${ASSEMBLY_PATH} does not exist. Skipping ${BARCODE_NAME}..."
    fi
done

echo "Minimap2 processing completed for all samples."
