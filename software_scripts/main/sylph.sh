#!/bin/bash

# Activar el entorno de conda para Sylph
echo "Activating conda environment for Sylph..."
# source ~/anaconda3/etc/profile.d/conda.sh
# conda activate sylph

# Definir rutas
FASTQ_DIR="/media/nova/datos/wastewaters_ch/results/combined_barcodes"
GTDB_DB="/media/nova/datos/databases/sylph/gtdb-r220-c200-dbv1.syldb"
# Salidas
SKETCH_DIR="/media/nova/datos/wastewaters_ch/results/sylph_results_combined_barcodes/sketch_dir"
OUTPUT_DIR="/media/nova/datos/wastewaters_ch/results/sylph_results_combined_barcodes/output_dir"

# # Crear directorios de salida si no existen
# mkdir -p "$OUTPUT_DIR"
# mkdir -p "$SKETCH_DIR"

# # Paso 1: Crear sketches para los archivos FASTQ
# echo "Sketching FASTQ files..."
# sylph sketch "$FASTQ_DIR"/*.fastq -d "$SKETCH_DIR"

# # Paso 2: Realizar el perfilado contra la base de datos GTDB-R220 sin opción '-u'
# echo "Profiling against GTDB database without '-u' option..."
# sylph profile "$GTDB_DB" "$SKETCH_DIR"/*.sylsp -t 10 -o "$OUTPUT_DIR/results_profile.tsv"

# # Paso 3: Realizar el perfilado con opción '-u'
# echo "Profiling against GTDB database with '-u' option..."
# sylph profile -u "$GTDB_DB" "$SKETCH_DIR"/*.sylsp -t 10 -o "$OUTPUT_DIR/results_with_unknown.tsv"

# # Paso 4: Realizar Query contra la base de datos GTDB-R220
# echo "Querying GTDB database..."
# sylph query "$GTDB_DB" "$SKETCH_DIR"/*.sylsp -t 10 -o "$OUTPUT_DIR/results_query.tsv"

# echo "Analysis complete. Results saved in $OUTPUT_DIR"

# Paso 5: Contar la cantidad de lecturas para cada archivo FASTQ y guardarlas en archivos TSV
echo "Counting reads for each FASTQ file..."

# Loop through each .fastq or .fastq.gz file
for FASTQ_FILE in "$FASTQ_DIR"/*.fastq*; do
    FILENAME=$(basename "$FASTQ_FILE" .fastq.gz)
    FILENAME=$(basename "$FILENAME" .fastq)  # Handle both .fastq and .fastq.gz
    echo "Counting for FASTQ file: $FASTQ_FILE"
    
    # Determine if the file is gzipped
    if [[ "$FASTQ_FILE" == *.gz ]]; then
        READ_COMMAND="zcat"  # Use zcat for .gz files
    else
        READ_COMMAND="cat"   # Use cat for uncompressed files
    fi

    # Verify that every fourth line starts with '@'
    FORMAT_VALID=1
    $READ_COMMAND "$FASTQ_FILE" | awk -v filename="$FILENAME" '
        NR % 4 == 1 {
            if ($0 !~ /^@/) {
                print filename "\tInvalid Format" >> "'"$OUTPUT_DIR/countreads.tsv"'"
                FORMAT_VALID=0
                exit 1
            }
        }
    '

    # Proceed with counting if the format is correct
    if [ "$FORMAT_VALID" -eq 1 ]; then
        READ_COUNT=$($READ_COMMAND "$FASTQ_FILE" | awk 'NR % 4 == 1' | wc -l)
        echo -e "$FILENAME\t$READ_COUNT" >> "$OUTPUT_DIR/countreads.tsv"
    else
        echo -e "$FILENAME\tInvalid Format" >> "$OUTPUT_DIR/countreads.tsv"
    fi
done

echo "Read counts saved in $OUTPUT_DIR/countreads.tsv"

#conda deactivate