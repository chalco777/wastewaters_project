#!/bin/bash

# Converts per-barcode BAMs from multiple pod5 basecalling runs into FASTQ,
# then merges same barcodes across runs.
# Input: base_path, expected to contain directories named like basecalling_pod5skip_*;
#   inside each, demultiplexed barcode subdirectories under results/demuxed/demuxed/
#   each with a reads.bam. 
# Output: per-pod per-barcode FASTQ files at:
#     <output_path>/basecalling_pod5skip_<...>/<barcode>.fastq.gz
#   and merged final FASTQ per barcode across pods at:
#     <output_path>/combined_fastq/final_<barcode>.fastq.gz

# Define absolute paths
base_path="/media/nova/datos/Diego/wastewater/240920_basecalling_p2_temporal"
output_path="/media/nova/datos/wastewaters_ch/results/fastq_temporal"
mkdir -p "$output_path"  # Create the main output directory if it doesn't exist

# Loop through each pod5skip directory
for pod_dir in "$base_path"/basecalling_pod5skip_*; do
  if [[ -d "$pod_dir" ]]; then
    pod_name=$(basename "$pod_dir")  # Get the pod5skip name
    pod_output_path="$output_path/$pod_name"  # Create a subfolder for this pod5skip
    mkdir -p "$pod_output_path"  # Create the output directory for this pod5skip if it doesn't exist
    demuxed_path="$pod_dir/results/demuxed/demuxed"

    # Loop through each barcode directory in the demuxed path
    for barcode_dir in "$demuxed_path"/*; do
      if [[ -d "$barcode_dir" ]]; then
        name=$(basename "$barcode_dir")
        barcode_name=$(echo "$name" | cut -d '_' -f3)
        bam_file="$barcode_dir/reads.bam"

        # Convert BAM to FASTQ using samtools
        if [[ -f "$bam_file" ]]; then
          fastq_output="$pod_output_path/${barcode_name}.fastq.gz"
          samtools fastq "$bam_file" | gzip > "$fastq_output"
          echo "Converted $bam_file to $fastq_output"
        else
          echo "No BAM file found in $barcode_dir"
        fi
      fi
    done
  fi
done

# Create a folder for the combined FASTQ files
combined_output="$output_path/combined_fastq"
mkdir -p "$combined_output"

# Combine FASTQ files of the same barcode from different pod5skip directories
for barcode_name in $(ls "$output_path"/basecalling_pod5skip_*/ | cut -d'.' -f1 | sort | uniq); do
  fastq_files=("$output_path"/basecalling_pod5skip_*/"${barcode_name}.fastq.gz")
  
  # Ensure there are FASTQ files to combine
  if [[ ${#fastq_files[@]} -gt 0 ]]; then
    combined_fastq="$combined_output/final_${barcode_name}.fastq.gz"
    cat "${fastq_files[@]}" > "$combined_fastq"
    echo "Combined FASTQ files for $barcode_name into $combined_fastq"
  else
    echo "No FASTQ files found for $barcode_name"
  fi
done