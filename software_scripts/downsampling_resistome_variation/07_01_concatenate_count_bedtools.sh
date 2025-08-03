#!/bin/bash
#Combines all coverage files (output from bedtools coverage) into a single file, 
#adding a column with the sample name for easy comparison then with the Rmd.
# Directorio con los archivos de coverage
COVERAGE_DIR="/media/crowfoot2/DATOS/wastewaters_ch/results/count_bedtools"

# Archivo de salida
OUTPUT_FILE="$COVERAGE_DIR/combined_coverage.tsv"

# Crear el archivo de salida y escribir el encabezado
echo -e "Muestra\tContig\tStart\tEnd\tFeature\tConteo\tNumero_bases_cubiertas\tLongitud_del_gen\tFraccion_cubierta" > "$OUTPUT_FILE"

# Iterar sobre cada archivo de coverage
for coverage_file in "$COVERAGE_DIR"/*_coverage.txt; do
    # Extraer el nombre de la muestra del archivo
    sample_name=$(basename "$coverage_file" | sed 's/_coverage\.txt//')

    # Agregar la columna de la muestra y concatenar al archivo de salida
    awk -v sample="$sample_name" '{print sample "\t" $0}' "$coverage_file" >> "$OUTPUT_FILE"
    echo "Procesado: $coverage_file"
done

echo "Archivo combinado creado en: $OUTPUT_FILE"
