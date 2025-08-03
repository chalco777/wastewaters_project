#!/bin/bash

## Script functionality:

# This script unifies scattered reads across three different dirs by barcode,
# concatenating existing .fastq.gz files to create a single combined
# barcode file in the output directory results/combined_barcodes. This is useful for consolidating
# fragmented sequencing data by tag prior to downstream analysis.


# Definir directorios de origen
dir1="/media/nova/datos/wastewaters_ch/results/fastq_temporal/filtlong"
dir2="/media/nova/datos/Diego/wastewater/240924_MA_p2_temporal-kracken-bracken/data/170524_HMA/fastq_pass"
dir3="/media/nova/datos/Diego/wastewater/240924_MA_p2_temporal-kracken-bracken/data/150524_HMA/fastq_pass"

# Directorio de salida
output_dir="/media/nova/datos/wastewaters_ch/results/combined_barcodes"

# Crear el directorio de salida si no existe
mkdir -p "$output_dir"

# Obtener lista de barcodes de los tres directorios
barcodes=()

# Función para extraer barcodes de un directorio
get_barcodes() {
    local dir="$1"
    local pattern="$2"
    for file in "$dir"/$pattern; do
        if [ -e "$file" ]; then
            filename=$(basename "$file")
            # Extraer el barcode del nombre del archivo usando expresión regular
            if [[ "$filename" =~ (barcode[0-9]+) ]]; then
                barcode="${BASH_REMATCH[1]}"
                # Añadir barcode a la lista si no está ya incluido
                if [[ ! " ${barcodes[@]} " =~ " ${barcode} " ]]; then
                    barcodes+=("$barcode")
                fi
            fi
        fi
    done
}

# Extraer barcodes de dir1
get_barcodes "$dir1" "final_barcode*.fastq.gz"

# Extraer barcodes de dir2
get_barcodes "$dir2" "barcode*_combined.fastq.gz"

# Extraer barcodes de dir3
get_barcodes "$dir3" "barcode*_combined.fastq.gz"

# Mostrar lista de barcodes encontrados
echo "Barcodes encontrados: ${barcodes[@]}"

# Iterar sobre cada barcode
for barcode in "${barcodes[@]}"; do
    echo "Procesando $barcode"

    # Inicializar una lista de archivos a concatenar
    files_to_concat=()

    # Buscar el archivo en dir1
    file1="$dir1/final_${barcode}.fastq.gz"
    if [ -f "$file1" ]; then
        files_to_concat+=("$file1")
    else
        echo "Advertencia: Archivo $file1 no encontrado."
    fi

    # Buscar el archivo en dir2
    file2="$dir2/${barcode}_combined.fastq.gz"
    if [ -f "$file2" ]; then
        files_to_concat+=("$file2")
    else
        echo "Advertencia: Archivo $file2 no encontrado."
    fi

    # Buscar el archivo en dir3
    file3="$dir3/${barcode}_combined.fastq.gz"
    if [ -f "$file3" ]; then
        files_to_concat+=("$file3")
    else
        echo "Advertencia: Archivo $file3 no encontrado."
    fi

    # Verificar si hay archivos para concatenar
    if [ ${#files_to_concat[@]} -gt 0 ]; then
        # Definir el archivo de salida
        output_file="$output_dir/${barcode}_final_combined.fastq.gz"

        # Concatenar los archivos
        cat "${files_to_concat[@]}" > "$output_file"

        echo "Archivos combinados para $barcode en $output_file"
    else
        echo "No se encontraron archivos para concatenar para $barcode."
    fi
done
