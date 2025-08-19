#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
Pipeline to classify combined per-barcode FASTQ reads with Kraken2 and
estimate abundances with Bracken.

Inputs:
  - params.fastqDir: directory containing combined barcode FASTQ files matching
      "*_final_combined.fastq.gz". In this case results/combined_barcodes.
  - params.kraken2DB: path to the Kraken2 database to use for classification.

Outputs:
  - Kraken2: per-sample classification and report files written to the directory
      published by the Kraken2 process (here: /media/nova/datos/wastewaters_ch/results/kraken2)
      as <sample_id>_kraken2_classification.tsv and <sample_id>_kraken2_report.tsv.
  - Bracken: per-sample abundance estimates and summary reports written to
      /media/nova/datos/wastewaters_ch/results/bracken as <sample_id>.bracken
      and <sample_id>.outreport_bracken.
*/


params.kraken2DB = "/media/nova/datos/k2_standard_20240605" 
params.fastqDir = "/media/nova/datos/wastewaters_ch/results/combined_barcodes"
params.kraken="/media/nova/datos/wastewaters_ch/results/kraken2"
process Kraken2 {
    conda '/home/nova/anaconda3/envs/kraken2'
    publishDir '/media/nova/datos/wastewaters_ch/results/kraken2', mode: 'copy' 
    cpus 15
    maxForks 1

    input:
    tuple val(sample_id), path(query_file)

    output:
    val sample_id
    path "${sample_id}_kraken2_report.tsv"
    path "${sample_id}_kraken2_classification.tsv"

    script:
    """
    kraken2 --db $params.kraken2DB \\
    --threads ${task.cpus} \\
    --gzip-compressed \\
    --report ${sample_id}_kraken2_report.tsv \\
    --output ${sample_id}_kraken2_classification.tsv \\
    $query_file
    """
}

process Bracken {
    conda '/home/nova/anaconda3/envs/kraken2'
    publishDir '/media/nova/datos/wastewaters_ch/results/bracken', mode: 'copy' 
    cpus 15
    maxForks 1

    input:
    val sample_id
    path query_file
    path classification

    output:
    val sample_id
    path "${sample_id}.outreport_bracken"
    path "${sample_id}.bracken"

    script:
    """
    bracken \\
    -d $params.kraken2DB \\
    -i $query_file \\
    -o ${sample_id}.bracken \\
    -w ${sample_id}.outreport_bracken \\
    -r 1000 \\
    -l S
    """

}

workflow {
    channel_fastq = channel.fromFilePairs("$params.fastqDir/*_final_combined.fastq.gz", size:1)
    channel_fastq.view()
    Kraken2(channel_fastq)
    Bracken(Kraken2.out)
}




