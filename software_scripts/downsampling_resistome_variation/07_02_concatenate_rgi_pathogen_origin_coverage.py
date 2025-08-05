#!/usr/bin/env python3

"""
Merges BED-pathogen-origin annotations with coverage results per barcode and
produce both per-barcode merged tables and a combined summary.

Inputs:
  - bed_patho_dir (results/bed_PATHO_results): directory with BED files annotated with pathogen-origin predictions,
      named like <barcode>.bed (columns: contig, start, end, annotation, species).
  - count_bedtools_dir (results/count_bedtools): directory with coverage outputs from bedtools, named like
      <barcode>_coverage.txt (columns from coverage: chromosome, start, end, feature,
      count, bases_covered, feature_length, fraction_covered).

OUTPUT DIR: /results/merged_countrgi_withpatho
Per-barcode output:
  - <barcode>_merged.txt in output_dir: inner join of coverage and BED annotations on
    (Chromosome, Start, End), retaining coverage metrics and pathogen annotation.

Combined output:
  - all_merged_combined.txt: all per-barcode merged tables concatenated with a leading
    "Muestra" column indicating the barcode/sample.

"""


import os
import pandas as pd
from glob import glob

# Directorios de entrada
bed_patho_dir = "/media/crowfoot2/DATOS/wastewaters_ch/results/bed_PATHO_results"
count_bedtools_dir = "/media/crowfoot2/DATOS/wastewaters_ch/results/count_bedtools"

# Directorio de salida
output_dir = "/media/crowfoot2/DATOS/wastewaters_ch/results/merged_countrgi_withpatho"
os.makedirs(output_dir, exist_ok=True)

# Obtener lista de archivos BED y de cobertura
bed_files = glob(os.path.join(bed_patho_dir, "*.bed"))
count_files = glob(os.path.join(count_bedtools_dir, "*_coverage.txt"))

# Crear diccionarios para mapear los archivos por código de barras
bed_files_dict = {}
for f in bed_files:
    basename = os.path.basename(f)
    barcode = basename.replace(".bed", "")
    bed_files_dict[barcode] = f

count_files_dict = {}
for f in count_files:
    basename = os.path.basename(f)
    barcode = basename.replace("_coverage.txt", "")
    count_files_dict[barcode] = f

# Procesar cada código de barras y guardar archivos intermedios
merged_files = []
for barcode in bed_files_dict:
    if barcode in count_files_dict:
        bed_file = bed_files_dict[barcode]
        count_file = count_files_dict[barcode]

        # Leer el archivo BED en un DataFrame
        bed_df = pd.read_csv(
            bed_file, sep='\t', header=None,
            names=['Chromosome', 'Start', 'End', 'Annotation', 'Pathogen']
        )

        # Leer el archivo de cobertura en un DataFrame
        count_df = pd.read_csv(
            count_file, sep='\t', header=None,
            names=['Chromosome', 'Start', 'End', 'Feature', 'Conteo',
                   'Numero_bases_cubiertas', 'Longitud_del_gen', 'Fraccion_cubierta']
        )

        # Convertir las columnas clave a tipos adecuados
        bed_df['Chromosome'] = bed_df['Chromosome'].astype(str)
        bed_df['Start'] = bed_df['Start'].astype(int)
        bed_df['End'] = bed_df['End'].astype(int)

        count_df['Chromosome'] = count_df['Chromosome'].astype(str)
        count_df['Start'] = count_df['Start'].astype(int)
        count_df['End'] = count_df['End'].astype(int)

        # Realizar la unión de los DataFrames basándose en 'Chromosome', 'Start' y 'End'
        merged_df = pd.merge(count_df, bed_df, on=['Chromosome', 'Start', 'End'], how='left')

        # Seleccionar y reordenar las columnas según sea necesario
        merged_df = merged_df[['Chromosome', 'Start', 'End', 'Feature', 'Conteo',
                               'Numero_bases_cubiertas', 'Longitud_del_gen',
                               'Fraccion_cubierta', 'Pathogen']]

        # Guardar en un archivo
        output_file = os.path.join(output_dir, f"{barcode}_merged.txt")
        merged_df.to_csv(output_file, sep='\t', index=False)
        merged_files.append(output_file)
        print(f"Archivo procesado y guardado: {output_file}")
    else:
        print(f"No se encontró archivo de cobertura para {barcode}")

# Combinar todos los archivos 'merged' en un solo DataFrame
all_merged_df = pd.DataFrame()
for file in merged_files:
    # Extraer el nombre de la muestra del nombre del archivo
    sample_name = os.path.basename(file).replace('_merged.txt', '')

    # Leer el archivo y agregar la columna 'Muestra'
    df = pd.read_csv(file, sep='\t')
    df.insert(0, 'Muestra', sample_name)  # Insertar la columna al inicio

    # Concatenar el DataFrame al acumulador
    all_merged_df = pd.concat([all_merged_df, df], ignore_index=True)

# Guardar el archivo final concatenado
final_output_file = os.path.join(output_dir, "all_merged_combined.txt")
all_merged_df.to_csv(final_output_file, sep='\t', index=False)
print(f"Archivo combinado final guardado: {final_output_file}")
