#!/bin/bash

# Runs Flye assemblies on downsampled FASTQ files (assumed Nanopore corrected reads)
# Input directory: contains per-barcode downsampled .fastq.gz files (BASE_DIR).
#   e.g. /media/crowfoot2/DATOS/wastewaters_ch/results/downsampling_result
# Output directory: MAIN_OUTPUT_DIR, with one subdirectory per barcode holding Flye results.
#   e.g. /media/crowfoot2/DATOS/wastewaters_ch/results/DS_flye/<barcode_name>/

# Crear el directorio principal de salida
MAIN_OUTPUT_DIR="/media/crowfoot2/DATOS/wastewaters_ch/results/DS_flye"
mkdir -p $MAIN_OUTPUT_DIR

# Directorio base donde están ubicados los archivos FASTQ comprimidos
BASE_DIR="/media/crowfoot2/DATOS/wastewaters_ch/results/downsampling_result"

# Activar el entorno de conda con Flye
source $(conda info --base)/etc/profile.d/conda.sh
conda activate flye

# Iterar sobre cada archivo FASTQ en el directorio base
for FASTQ_FILE in $BASE_DIR/*.fastq.gz; do
    # Obtener el nombre del archivo sin la extensión (.fastq.gz)
    BARCODE_NAME=$(basename "$FASTQ_FILE" .fastq.gz)
    
    # Crear el directorio de salida para el ensamblaje de este archivo
    BARCODE_OUTPUT_DIR="${MAIN_OUTPUT_DIR}/${BARCODE_NAME}"
    mkdir -p "$BARCODE_OUTPUT_DIR"

    # Ejecutar Flye en el archivo FASTQ actual
    flye --nano-corr "$FASTQ_FILE" -o "$BARCODE_OUTPUT_DIR" -t 30 --meta

    echo "Ensamblaje completado para $BARCODE_NAME"
done

echo "Flye assembly completed for all FASTQ samples."
