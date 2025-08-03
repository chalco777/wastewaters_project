#!/bin/bash


## Script functionality:

# Filters compressed FASTQ files from a specified input directory using filtlong,
# retaining reads of at least $min_length (200 bases in this case). 
# Input directory: contains combined .fastq.gz files to be filtered (here: 
#   /media/nova/datos/wastewaters_ch/results/fastq_temporal/combined_fastq).
# Output directory: receives the filtered, recompressed FASTQ files with the same
# base names (here: /media/nova/datos/wastewaters_ch/results/fastq_temporal/filtlong).



# Definir rutas absolutas
input_dir="/media/nova/datos/wastewaters_ch/results/fastq_temporal/combined_fastq"
output_dir="/media/nova/datos/wastewaters_ch/results/fastq_temporal/filtlong"
min_length=200

# Crear el directorio de salida si no existe
mkdir -p "$output_dir"

# Bucle a través de todos los archivos FASTQ en el directorio de entrada
for fastq_file in "$input_dir"/*.fastq.gz; do
  # Obtener el nombre base del archivo
  base_name=$(basename "$fastq_file")

  # Definir el nombre del archivo de salida
  output_file="$output_dir/$base_name"

  # Aplicar Filtlong con el largo mínimo
  filtlong --min_length $min_length "$fastq_file" | gzip > "$output_file"

  echo "Filtrado de $fastq_file completado. Resultado guardado en $output_file"
done
