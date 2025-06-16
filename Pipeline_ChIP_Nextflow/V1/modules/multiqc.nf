#!/usr/bin/env nextflow


/*
    Launch MultiQC report for a folder including FastQC reports
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
