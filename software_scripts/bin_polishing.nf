#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.fastqDir = "$launchDir/reads"
params.assembly= "$launchDir/preliminarbins"

process Mapwithminimap2 {
    conda '/home/crowfoot2/anaconda3/envs/medaka'
    cpus=2
    publishDir 'results/minimap2', mode: 'copy'
    input:
    tuple val(sample_id), path(query_file)
    path assembly
    each number
    output:
    val sample_id
    path "${sample_id}_${number}.bam"
    val number
    script:
    """
    minimap2 -ax map-ont $assembly/bin.${number}.fa $query_file \\
        | samtools sort \\
        | samtools view -F 2308 -b > ${sample_id}_${number}.bam
    """
}
process FromBAM_toFastq {
    conda '/home/crowfoot2/anaconda3/envs/medaka'
    publishDir 'results/fastq', mode: 'copy'
    input:
    val sample_id
    path query_file
    val number
    output:
    val sample_id
    path "${sample_id}_${number}.fastq"
    val number
    script:
    """
    samtools fastq $query_file > ${sample_id}_${number}.fastq
    """
}
process Polishing_medaka {
    conda '/home/crowfoot2/anaconda3/envs/medaka_gpu'
    cpus=30
    publishDir 'results/medaka_gpu', mode: 'copy'
    input:
    val sample_id
    path read
    val number
    path assembly
    output:
    val sample_id
    path "medaka.${sample_id}.${number}"
    val number
    script:
    """
    medaka_consensus -i $read \\
    -d $assembly/bin.${number}.fa \\
    -o medaka.${sample_id}.${number} \\
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
    Mapwithminimap2(channel_fastq,params.assembly,number)
    FromBAM_toFastq(Mapwithminimap2.out)
    Polishing_medaka(FromBAM_toFastq.out,params.assembly)
    Checkm2_polished(Polishing_medaka.out)
}
