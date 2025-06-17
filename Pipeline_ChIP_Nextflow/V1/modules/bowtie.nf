#!/usr/bin/env nextflow


/*
    Run bowtie2 on the trim-illumina.fastq.gz
*/

process bowtieIndex {
    publishDir 'data', mode: 'copy'

    input:
        path genomeFasta
        val batch

    output:
        path "${batch}/genome_index.*.bt2"
        path "${batch}", emit: indexDir

    script:
    """
    mkdir -p ${batch}
    bowtie2-build $genomeFasta ${batch}/genome_index
    """
}


process bowtieAlign {
    publishDir 'results', mode: 'copy'
    
    input: 
        path trimDir
        path genomeIndex
        val threads
        val batch

    output:
        path "${batch}_bowtieAlign/*", emit : bowtieDir
        path "${batch}_bowtieAlign/*/*sam"

    script:
    """
    name=\$(basename ${trimDir})
    r1=\$(find -L "${trimDir}" -type f -name "*R1*.fastq*" | sort | head -n 1)
    r2=\$(find -L "${trimDir}" -type f -name "*R2*.fastq*" | sort | head -n 1)
    mkdir -p ${batch}_bowtieAlign/\${name}
    bowtie2 \\
        -p $threads  \\
        --very-sensitive \\
        --phred33 \\
        --no-mixed \\
        --no-discordant \\
        --dovetail \\
        -x ${genomeIndex}/genome_index \\
        -1 \$r1 -2 \$r2 \\
        -S ${batch}_bowtieAlign/\${name}/\${name}.sam
    """
}