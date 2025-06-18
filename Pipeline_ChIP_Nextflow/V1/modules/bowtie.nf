#!/usr/bin/env nextflow


/*
    Run bowtie2 on the trim-illumina.fastq
*/

process bowtieIndex {
    publishDir 'data', mode: 'copy'

    input:
        path genomeFasta
        val batch

    output:
        path "${batch}", emit: bowtieIndexDir

    script:
    """
    mkdir -p ${batch}
    bowtie2-build $genomeFasta ${batch}/genome_index
    """
}


process bowtieAlign {
    input: 
        tuple val (name), path (r1), path (r2)
        path bowtieIndexDir
        val threads
        val batch

    output:
        tuple val (name), path ("${batch}_4-bowtieAlign/${name}/${name}.sam"), emit : bowtieAligned

    script:
    """
    mkdir -p ${batch}_4-bowtieAlign/${name}
    bowtie2 \\
        -p $threads  \\
        --very-sensitive \\
        --phred33 \\
        --no-mixed \\
        --no-discordant \\
        --dovetail \\
        -x ${bowtieIndexDir}/genome_index \\
        -1 $r1 -2 $r2 \\
        -S ${batch}_4-bowtieAlign/${name}/${name}.sam
    """
}