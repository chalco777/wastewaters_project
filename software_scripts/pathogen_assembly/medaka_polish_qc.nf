#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.fastqDir = "$launchDir/reads"
params.assembly= "$launchDir/preliminarbins"
process Polishing_medaka {
    conda '/home/crowfoot2/anaconda3/envs/medaka_gpu'
    cpus=30
    publishDir 'results/medaka', mode: 'copy'
    input:
    tuple val(sample_id), path(query_file)
    path assembly
    each number
    output:
    val sample_id
    path "medaka.${sample_id}_${number}"
    val number
    script:
    """
    medaka_consensus -i $query_file \\
    -d $assembly/bin.${number}.fa \\
    -o medaka.${sample_id}_${number} \\
    -t ${task.cpus} -b 50
    """
}
process Checkm2_polished {
    conda '/home/crowfoot2/anaconda3/envs/checkm2'
    cpus=10
    publishDir 'results/checkm2_polished', mode: 'copy'
    input:
    val sample_id
    path query_file
    val number
    output:
    val sample_id
    path "polished_${sample_id}_${number}"
    script:
    """
    checkm2 predict --threads ${task.cpus} -x .fasta\\
     --input ${query_file}/consensus.fasta\\
     --output-directory polished_${sample_id}_${number}
    """
}
workflow {
    channel_fastq = channel.fromFilePairs("${params.fastqDir}/*.fastq", size:1)
    number=[1,2,3,4,5,6,7]
    channel_fastq.view()
    Polishing_medaka(channel_fastq, params.assembly,number)
    Checkm2_polished(Polishing_medaka.out)
}
