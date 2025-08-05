#!/bin/bash

# Computes coverage of RGI-predicted regions (BED) over aligned reads for each barcode.
# Input:
#   - BED_DIR (results/rgi_D_beds): BED files with resistance gene regions (e.g., from RGI-to-BED conversion), expected to match pattern "*GB.bed".
#   - ALIGNMENT_DIR (results/minimap_result): per-barcode sorted BAMs from minimap2 alignment (named like <barcode>_sorted.bam under <barcode>/).
#   - FAI source (results/DS_flye): Flye assemblies in FAI format under FAI_DIR/<barcode>/assembly.fasta[.fai].
# Output:
#   - OUTPUT_DIR (results/count_bedtools): for each barcode, produces:
#       * <barcode>_genome.txt : genome file (contig name and length) derived from the FASTA index.
#       * <barcode>_coverage.txt : coverage of each BED interval against the BAM, using bedtools coverage.
# The script ensures required indices exist, sorts BEDs if needed, and applies a 50% fraction threshold (-f 0.5).


# Directorios base
BED_DIR="/media/crowfoot2/DATOS/wastewaters_ch/results/rgi_D_beds"
ALIGNMENT_DIR="/media/crowfoot2/DATOS/wastewaters_ch/results/minimap_result"
FAI_DIR="/media/crowfoot2/DATOS/wastewaters_ch/results/DS_flye"
OUTPUT_DIR="/media/crowfoot2/DATOS/wastewaters_ch/results/count_bedtools"

# Crear el directorio de salida si no existe
mkdir -p "$OUTPUT_DIR"

# Iterar sobre cada archivo .bed
for bed_file in "$BED_DIR"/*GB.bed; do
    # Obtener el nombre del archivo sin extensión
    filename=$(basename "$bed_file")
    base_name="${filename%.bed}"  # Eliminar la extensión .bed
    sorted_bed_file="$BED_DIR/${base_name}.sorted.bed"
    sort -k1,1 -k2,2n "$bed_file" > "$sorted_bed_file"
    echo "Archivo BED ordenado: $sorted_bed_file"
    # Buscar el directorio del archivo BAM correspondiente
    # Suponemos que el nombre del directorio es similar al base_name pero reemplazando el sufijo de tamaño por 1405Gb
    bam_dir="$ALIGNMENT_DIR/$base_name"
    bam_file="$bam_dir/${base_name}_sorted.bam"

    # Verificar si el archivo BAM existe antes de ejecutar bedtools
    if [[ -f "$bam_file" ]]; then
        # Generar el genome file usando el archivo assembly.fasta.fai correspondiente
        fai_dir="$FAI_DIR/$base_name"
        fai_file="$fai_dir/assembly.fasta.fai"
        genome_file="$OUTPUT_DIR/${base_name}_genome.txt"

        # Si el archivo .fai no existe, crearlo con samtools faidx
        if [[ ! -f "$fai_file" ]]; then
            fasta_file="$fai_dir/assembly.fasta"
            if [[ -f "$fasta_file" ]]; then
                samtools faidx "$fasta_file"
                echo "Archivo FAI creado: $fai_file"
            else
                echo "Archivo FASTA no encontrado para $bam_base_name en $fai_dir. Saltando..."
                continue
            fi
        fi

        # Crear el genome file si el FAI existe
        awk '{print $1"\t"$2}' "$fai_file" > "$genome_file"
        echo "Genome file creado: $genome_file"

        # Ejecutar bedtools coverage y guardar el resultado en la carpeta de salida
        output_file="$OUTPUT_DIR/${base_name}_coverage.txt"
        bedtools coverage -a "$sorted_bed_file" -b "$bam_file" -g "$genome_file" -f 0.5 -sorted > "$output_file"
        
        echo "Procesado: $sorted_bed_file con $bam_file y $genome_file -> $output_file"
    else
        echo "Archivo BAM no encontrado para $base_name en $bam_dir. Saltando..."
    fi

    echo "El fin"
done
