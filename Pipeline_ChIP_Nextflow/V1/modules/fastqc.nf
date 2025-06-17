#!/usr/bin/env nextflow


/*
    Launch FastQC report for a tuple of R1, R2
*/

process reportRawFastqc {
    publishDir 'results', mode: 'copy'
    
    input:
        tuple val (name), path (r1), path (r2)
        val batch

    output:
        path "${batch}_rawReports/fastQC_${name}", emit : fastqcDir
        path "${batch}_rawReports/fastQC_${name}/*.zip"
        path "${batch}_rawReports/fastQC_${name}/*.html"

    script:
    """
    mkdir -p ${batch}_rawReports/fastQC_${name}
    fastqc ${r1} ${r2} -o ${batch}_rawReports/fastQC_${name}
    """
}

process reportTrimFastqc {
    publishDir 'results', mode: 'copy'
    
    input:
        path trimDir
        val batch

    output:
        path "${batch}_trimReports/fastQC_*", emit : fastqcDir
        path "${batch}_trimReports/fastQC_*/*.zip"
        path "${batch}_trimReports/fastQC_*/*.html"

    script:
    """
    name=\$(basename ${trimDir})
    r1=\$(find -L "${trimDir}" -type f -name "*R1*.fastq*" | sort | head -n 1)
    r2=\$(find -L "${trimDir}" -type f -name "*R2*.fastq*" | sort | head -n 1)
    mkdir -p ${batch}_trimReports/fastQC_\${name}
    fastqc \${r1} \${r2} -o ${batch}_trimReports/fastQC_\${name}
    """
}