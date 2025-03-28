#!/usr/bin/env nextflow

/*
 * Create a MultiQC report for the PE reads contained 
 * in the text file given as arguments.
 */

nextflow.enable.dsl=2

// Display help message
if (params.help) {
    log.info """
    Usage: nextflow run fastqc.nf --fastq_list <path> [--output_fastqc <path>] [--output_multiqc <path>] [--help]

    Required parameters:
    --fastq_list    Path to a text file containing paired-end FASTQ files to be quality checked

    Optional parameters:
    --output_fastqc  Output directory for the FastQC reports (default: ./01_FastQC-report/data)
    --output_multiqc Output directory for the MultiQC report (default: ./01_FastQC-report)
    --help           Show this help message
    """
    exit 0
}

// Validate required parameters
if (!params.fastq_list) {
    exit 1, "ERROR: --fastq_list must be specified"
}

// Pipeline default parameters
params.output_fastqc = params.get('output_fastqc', "${params.output_root}/01_FastQC-report/data")           // FastQC Output directory
params.output_multiqc = params.get('output_multiqc', "${params.output_root}/01_FastQC-report")      // MultiQC Output directory

// FastQC reports, quality check of the raw fastq files
process QUALITY_CHECK {
    tag "FastQC of raw paired-end FASTQ files"

    input: 
        path reads_list from file(params.fastq_list)    // Illumina PE fastq list

    output:
        path("${params.output_fastqc}") into fastqc_reports
    
    script: 
    """
    mkdir -p ${params.output_fastqc}
    while read reads1 reads2; do
        fastqc -o ${params.output_fastqc} -t ${params.threads} ${reads1} ${reads2}
    done < <(cut -f1,2 -d$'\\t' $reads_list) 
    """
}

// MultiQC report, it merges FastQC reports into only one report
process GROUP_QC {
    tag "Combining FastQC reports"
    
    input:
        path fastqc_reports_dir from fastqc_reports.collect()         // FastQC reports
    
    output:
        path("${params.output_multiqc}/multiqc_data")
        path("${params.output_multiqc}/multiqc_report.html")

    script: 
    """
    multiqc ${fastqc_reports_dir} -o "${params.output_multiqc}"
    """

}

workflow {
    // Launch FastQC
    fastqc_reports = QUALITY_CHECK()
    // Launch MultiQC
    fastqc_reports | GROUP_QC()
}