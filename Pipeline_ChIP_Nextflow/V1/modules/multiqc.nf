#!/usr/bin/env nextflow


/*
    Launch MultiQC report for a folder including FastQC reports
    conda create -n multiqc_env python=3.10 multiqc=1.14 numpy=1.23
*/

process reportRawMultiqc {
    
    input:
        path fastqcDirs
        val batch

    output:
        path "${batch}_1-rawReports/multiqc_report.html"

    script:
    """
    conda run -n multiqc_env multiqc ${fastqcDirs} -o ${batch}_1-rawReports
    """
}

process reportTrimMultiqc {
    
    input:
        path fastqcDirs
        val batch

    output:
        path "${batch}_3-trimReports/multiqc_report.html"

    script:
    """
    conda run -n multiqc_env multiqc ${fastqcDirs} -o ${batch}_3-trimReports
    """
}
