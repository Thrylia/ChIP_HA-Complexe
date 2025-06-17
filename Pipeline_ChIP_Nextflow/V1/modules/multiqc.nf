#!/usr/bin/env nextflow


/*
    Launch MultiQC report for a folder including FastQC reports
    conda create -n multiqc_env python=3.10 multiqc=1.14 numpy=1.23
*/

process reportRawMultiqc {
    publishDir 'results', mode: 'copy'
    
    input:
        path fastqcDirs
        val batch

    output:
        path "${batch}_rawReports/multiqc_report.html"

    script:
    """
    conda run -n multiqc_env multiqc ${fastqcDirs} -o ${batch}_rawReports
    """
}

process reportTrimMultiqc {
    publishDir 'results', mode: 'copy'
    
    input:
        path fastqcDirs
        val batch

    output:
        path "${batch}_trimReports/multiqc_report.html"

    script:
    """
    conda run -n multiqc_env multiqc ${fastqcDirs} -o ${batch}_trimReports
    """
}
