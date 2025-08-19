#!/bin/bash

# Activate the conda environment for Sylph. We do this so Sylph commands are available.
echo "Activating conda environment for Sylph..."
source ~/anaconda3/etc/profile.d/conda.sh
conda activate sylph

# Set the paths for input and sylph gtdb database.
FASTQ_DIR="/media/nova/datos/wastewaters_ch/results/combined_barcodes"
GTDB_DB="/media/nova/datos/databases/sylph/gtdb-r220-c200-dbv1.syldb"

# Define output directories
SKETCH_DIR="/media/nova/datos/wastewaters_ch/results/sylph_results_combined_barcodes/sketch_dir"
OUTPUT_DIR="/media/nova/datos/wastewaters_ch/results/sylph_results_combined_barcodes/output_dir"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$SKETCH_DIR"

### SYLPH ANALYSIS STEPS ###

# Step 1: Generate sketches for the FASTQ files. We do this to prepare files for profiling.
echo "Sketching FASTQ files..."
sylph sketch "$FASTQ_DIR"/*.fastq -d "$SKETCH_DIR"

# Step 2: Profile against the GTDB-R220 database without the '-u' option.
echo "Profiling against GTDB database without '-u' option..."
sylph profile "$GTDB_DB" "$SKETCH_DIR"/*.sylsp -t 10 -o "$OUTPUT_DIR/results_profile.tsv"

# Step 3: Profile with the '-u' option enabled.
echo "Profiling against GTDB database with '-u' option..."
sylph profile -u "$GTDB_DB" "$SKETCH_DIR"/*.sylsp -t 10 -o "$OUTPUT_DIR/results_with_unknown.tsv"

# Step 4: Run a query against the GTDB-R220 database.
echo "Querying GTDB database..."
sylph query "$GTDB_DB" "$SKETCH_DIR"/*.sylsp -t 10 -o "$OUTPUT_DIR/results_query.tsv"

echo "Analysis complete. Results saved in $OUTPUT_DIR"

# Step 5: Count the number of reads for each FASTQ file and save the results in TSV files.
echo "Counting reads for each FASTQ file..."

# Loop through each .fastq file.
for FASTQ_FILE in "$FASTQ_DIR"/*.fastq*; do
    FILENAME=$(basename "$FASTQ_FILE" .fastq.gz)
    FILENAME=$(basename "$FILENAME" .fastq)  # Handle both .fastq and .fastq.gz
    echo "Counting for FASTQ file: $FASTQ_FILE"
    
    # Check if the file is gzipped.
    if [[ "$FASTQ_FILE" == *.gz ]]; then
        READ_COMMAND="zcat"  # Here use zcat for compressed file
    else
        READ_COMMAND="cat"   # And here use cat fr uncompressed files.
    fi

    # Make sure every fourth line starts with '@', we are doing this to validate the format
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

    # If the format is correct, proceed with counting.
    if [ "$FORMAT_VALID" -eq 1 ]; then
        READ_COUNT=$($READ_COMMAND "$FASTQ_FILE" | awk 'NR % 4 == 1' | wc -l)
        echo -e "$FILENAME\t$READ_COUNT" >> "$OUTPUT_DIR/countreads.tsv"
    else
        echo -e "$FILENAME\tInvalid Format" >> "$OUTPUT_DIR/countreads.tsv"
    fi
done

echo "Read counts saved in $OUTPUT_DIR/countreads.tsv"

conda deactivate