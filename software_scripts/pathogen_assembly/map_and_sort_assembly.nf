#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.baseDir = "/home/tao/Documents/alen-belen/concatenated_fastq"
params.assembly = "/home/tao/Documents/alen-belen/medaka_result"
params.outputDir = "/home/tao/Documents/alen-belen/minimap_result"
process Map_to_assembly {
    //si uso un proceso con conda todos conda, si uso docker
    //conda '/home/tao/anaconda3/envs/medaka'  // Especifica el entorno conda adecuado
    container "quay.io/biocontainers/minimap2:2.27--he4a0461_1"
    publishDir "${params.outputDir}", mode: 'copy'
    cpus 15
    maxForks 1

    input:
    tuple val(sample_id), path(query_file)
    path assembly

    output:
    path "${sample_id}_retroalignment.sam"
    val sample_id

    script:
    """
    minimap2 -ax map-ont -t $task.cpus ${assembly}/${sample_id}/consensus.fasta ${query_file} > ${sample_id}_retroalignment.sam
    """
}

process Sort_Samtools {
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
    publishDir "${params.outputDir}", mode: 'copy'
    cpus 14
    maxForks 1

    input:
    path sam
    val sample_id
    
    output:
    path "${sample_id}_retroalignment.bam"
    path "${sample_id}_retroalignment.bam.bai"
    path "${sample_id}__retroalignment_flagstat.txt"

    script:
    """
    samtools sort $sam \\
        --threads $task.cpus \\
        | samtools view -F 2048 -b > ${sample_id}_retroalignment.bam
    samtools index ${sample_id}_retroalignment.bam
    samtools flagstat ${sample_id}_retroalignment.bam > "${sample_id}_retroalignment_flagstat.txt"
    """
}

workflow {
    channel_fastq = channel.fromFilePairs("${params.baseDir}/barcode*", size:1)
    channel_fastq = channel.fromFilePairs("${params.baseDir}/barcode*", size:1)
    channel_fastq.view()
    Map_to_assembly(channel_fastq, params.assembly)
    Sort_Samtools(Map_to_assembly.out)

}


