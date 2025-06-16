#!/usr/bin/env nextflow


/*
    Launch FastQC report for a tuple of R1, R2
*/

process reportFastqc {
    publishDir 'results', mode: 'copy'
    
    input:
        tuple path (R1), path (R2)
        val batch

    output:
        path "${batch}_fastqcReport", emit : folderOut
        path "${batch}_fastqcReport/*.zip"
        path "${batch}_fastqcReport/*.html"

    script:
    """
    mkdir -p ${batch}_fastqcReport
    fastqc ${R1} ${R2} -o ${batch}_fastqcReport
    """
}
