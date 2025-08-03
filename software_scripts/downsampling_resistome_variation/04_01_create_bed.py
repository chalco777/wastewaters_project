#!/usr/bin/env python3

"""
Converts RGI "main" output tab-delimited summary (.txt) files into BED format.
Input:
  - BASE_DIR: directory (results/DS_rgi_result) containing per-sample RGI results organized as subdirectories
      named like "DS_<...>/"; inside each there are files like "barcode*.txt" produced
      by RGI main (with columns including Contig, Start, Stop, Best_Hit_ARO, Drug Class).
Output:
  - For each barcode .txt, writes a corresponding .bed file into OUTPUT_DIR (results/rgi_D_beds). Each BED line:
      contig (normalized), start (0-based), stop (end coordinate as in RGI), name combining
      Best_Hit_ARO and Drug Class with spaces replaced by underscores.
Usage example:
  ./rgi_to_bed.py \
    --base-dir /home/crowfoot2/DATOS/wastewaters_ch/results/DS_rgi_result \
    --output-dir /home/crowfoot2/DATOS/wastewaters_ch/results/rgi_D_beds
"""
import os
import csv

# Define el directorio base usando una ruta absoluta
base_dir = '/media/crowfoot2/DATOS/wastewaters_ch/results/DS_rgi_result'

# Define el directorio de salida para los archivos BED
output_dir = '/media/crowfoot2/DATOS/wastewaters_ch/results/rgi_D_beds'
os.makedirs(output_dir, exist_ok=True)

# Itera sobre cada directorio DS_* en el directorio base
for ds_dir in os.listdir(base_dir):
    ds_path = os.path.join(base_dir, ds_dir)
    if os.path.isdir(ds_path) and ds_dir.startswith('DS_'):
        print(f"Procesando directorio: {ds_path}")
        # Procesa cada archivo barcode*.txt en el directorio DS_*
        for file_name in os.listdir(ds_path):
            if file_name.endswith('.txt') and file_name.startswith('barcode'):
                barcode_file_path = os.path.join(ds_path, file_name)
                print(f"Procesando archivo: {barcode_file_path}")
                # Abre el archivo de código de barras
                with open(barcode_file_path, 'r') as infile:
                    reader = csv.reader(infile, delimiter='\t')
                    try:
                        header = next(reader)
                    except StopIteration:
                        print(f"Archivo vacío: {barcode_file_path}")
                        continue
                    # Obtiene los índices de los campos requeridos
                    try:
                        contig_idx = header.index('Contig')
                        start_idx = header.index('Start')
                        stop_idx = header.index('Stop')
                        best_hit_aro_idx = header.index('Best_Hit_ARO')
                        drug_class_idx = header.index('Drug Class')
                    except ValueError as e:
                        print(f"Error: Campo requerido faltante en el encabezado del archivo {barcode_file_path}: {e}")
                        continue
                    # Abre el archivo BED para escribir en el nuevo directorio de salida
                    bed_file_name = file_name.replace('.txt', '.bed')
                    bed_file_path = os.path.join(output_dir, bed_file_name)
                    with open(bed_file_path, 'w') as bedfile:
                        writer = csv.writer(bedfile, delimiter='\t', lineterminator='\n')
                        # Procesa cada línea de datos
                        for row in reader:
                            try:
                                contig = row[contig_idx]
                                # Elimina todo después e incluyendo el último guion bajo en el nombre del contig
                                if '_' in contig:
                                    contig = contig.rsplit('_', 1)[0]
                                start = int(row[start_idx]) - 1  # Ajusta para el formato BED (0-based)
                                stop = row[stop_idx]
                                best_hit_aro = row[best_hit_aro_idx]
                                drug_class = row[drug_class_idx].replace(' ', '_')                                # Combina Best_Hit_ARO y Drug Class con un guion bajo
                                name = f"{best_hit_aro}_{drug_class}"
                                writer.writerow([contig, start, stop, name])
                            except IndexError as e:
                                print(f"Error: Datos faltantes en el archivo {barcode_file_path}: {e}")
                            except ValueError as e:
                                print(f"Error: Datos inválidos en el archivo {barcode_file_path}: {e}")








