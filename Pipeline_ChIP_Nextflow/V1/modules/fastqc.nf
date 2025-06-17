#!/usr/bin/env nextflow


/*
    Launch FastQC report for a tuple of R1, R2
*/

process reportRawFastqc {
    
    input:
        tuple val (name), path (r1), path (r2)
        val batch

    output:
        path "${batch}_1-rawReports/fastQC_${name}", emit : fastqcDir
        path "${batch}_1-rawReports/fastQC_${name}/*.zip"
        path "${batch}_1-rawReports/fastQC_${name}/*.html"

    script:
    """
    mkdir -p ${batch}_1-rawReports/fastQC_${name}
    fastqc ${r1} ${r2} -o ${batch}_1-rawReports/fastQC_${name}
    """
}

process reportTrimFastqc {
    
    input:
        tuple val (name), path (r1), path (r2)
        val batch

    output:
        path "${batch}_3-trimReports/fastQC_*", emit : fastqcDir
        path "${batch}_3-trimReports/fastQC_*/*.zip"
        path "${batch}_3-trimReports/fastQC_*/*.html"

    script:
    """
    mkdir -p ${batch}_3-trimReports/fastQC_${name}
    fastqc ${r1} ${r2} -o ${batch}_3-trimReports/fastQC_${name}
    """
}