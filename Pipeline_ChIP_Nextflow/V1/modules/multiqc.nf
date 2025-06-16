#!/usr/bin/env nextflow


/*
    Launch MultiQC report for a folder including FastQC reports
    conda create -n multiqc_env python=3.10 multiqc=1.14 numpy=1.23
*/

process reportMultiqc {
    publishDir 'results', mode: 'copy'
    
    input:
        path fastqcReportFolder
        val batch

    output:
        path "${batch}_multiqcReport/multiqc_report.html"

    script:
    """
    conda run -n multiqc_env multiqc ${fastqcReportFolder} -o ${batch}_multiqcReport
    """
}
