#!/usr/bin/env python3

import os
import csv

"""
Combine RGI main output with RGI k-mer pathogen-of-origin predictions to produce
per-barcode BED files annotated with resistance gene info and predicted source species.

Inputs:
  - RGI_NORMAL_DIR (results/DS_rgi_result): contains per-DS subdirectories (e.g., DS_*) with RGI main summary
        .txt files (barcode*.txt) that include columns like ORF_ID, Contig, Start, Stop,
        Best_Hit_ARO, and Drug Class.
  - RGI_PATHO_DIR (results/rgi_PATHO_result): contains corresponding per-DS kmer_query results, with files named
        like <barcode>_61mer_analysis_rgi_summary.txt that include ORF_ID and a
        "CARD*kmer Prediction" field listing semicolon-separated predicted species.
Output (/results/bed_PATHO_results'):
  - For each barcode, writes a BED-style file to OUTPUT_DIR/<barcode>.bed with columns:
        contig (normalized), start (0-based), end, annotation (Best_Hit_ARO; Drug_Class),
        and species (each species produces its own line if multiple predicted per ORF).

Workflow:
  1. Parse pathogen-origin predictions per ORF from the k-mer summary.
  2. For each RGI main .txt, emit BED entries only for ORFs with pathogen predictions,
     duplicating lines when multiple species are associated.
"""

# Directorios base (actualiza estas rutas según tu sistema)
RGI_NORMAL_DIR = '/media/crowfoot2/DATOS/wastewaters_ch/results/DS_rgi_result'
RGI_PATHO_DIR = '/media/crowfoot2/DATOS/wastewaters_ch/results/rgi_PATHO_result'
OUTPUT_DIR = '/media/crowfoot2/DATOS/wastewaters_ch/results/bed_PATHO_results'

# Crear directorio de salida si no existe
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Función para obtener información de patógeno desde RGI Pathogen
def get_pathogen_info(ds_name, barcode_name):
    # Buscar el archivo correspondiente en RGI Pathogen
    ds_path = os.path.join(RGI_PATHO_DIR, ds_name)
    file_name = f"{barcode_name}_61mer_analysis_rgi_summary.txt"
    file_path = os.path.join(ds_path, file_name)

    pathogen_info = {}

    if os.path.exists(file_path):
        with open(file_path, 'r') as infile:
            reader = csv.DictReader(infile, delimiter='\t')
            for row in reader:
                orf_id = row['ORF_ID']
                card_kmer_prediction = row['CARD*kmer Prediction']
                # Obtener especies detectadas en 'CARD*kmer Prediction'
                species = card_kmer_prediction.split(';')
                species = [s.strip() for s in species]
                if orf_id not in pathogen_info:
                    pathogen_info[orf_id] = set()
                pathogen_info[orf_id].update(species)
    else:
        print(f"WARNING: No se encontró el archivo {file_path}")
    return pathogen_info

# Función principal
def main():
    for ds_name in os.listdir(RGI_NORMAL_DIR):
        ds_path = os.path.join(RGI_NORMAL_DIR, ds_name)
        if os.path.isdir(ds_path):
            for file_name in os.listdir(ds_path):
                if file_name.endswith('.txt'):
                    barcode_name = file_name.replace('.txt', '')
                    file_path = os.path.join(ds_path, file_name)

                    # Obtener información de patógeno
                    pathogen_info = get_pathogen_info(ds_name, barcode_name)

                    # Ruta del archivo BED por barcode
                    bed_file_path = os.path.join(OUTPUT_DIR, f'{barcode_name}.bed')

                    # Procesar el archivo .txt
                    with open(file_path, 'r') as infile:
                        reader = csv.DictReader(infile, delimiter='\t')
                        with open(bed_file_path, 'w') as outfile:
                            for row in reader:
                                orf_id = row['ORF_ID']
                                if orf_id in pathogen_info:
                                    contig = row['Contig']
                                    start = int(row['Start']) - 1  # Ajuste para el formato BED (restar 1)
                                    end = row['Stop']
                                    best_hit_aro = row['Best_Hit_ARO']
                                    drug_class = row['Drug Class']

                                    # Modificar el nombre del contig
                                    parts = contig.split('_')
                                    if len(parts) > 2:
                                        contig = '_'.join(parts[:2])

                                    # Información adicional para la cuarta columna del BED
                                    for species in pathogen_info[orf_id]:
                                        bed_info = f"{best_hit_aro}; {drug_class}"

                                        # Escribir las cuatro columnas en el archivo BED
                                        outfile.write(f"{contig}\t{start}\t{end}\t{bed_info}\t{species}\n")
    print("Proceso completado.")

if __name__ == "__main__":
    main()
