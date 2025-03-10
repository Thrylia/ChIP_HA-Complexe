#!/usr/bin/env nextflow

 /*
 * Marks and removes duplicates from a sorted bam 
 * file using Picard's Markduplicates package, and 
 * write a metrics file.
 */

nextflow.enable.dsl=2

// Display help message
if (params.help) {
    log.info """
    Usage: nextflow run script.nf --input_dir <path> [--output_dir <path>] [--help]

    Required parameters:
    --input_dir     Path to the directory containing sorted BAM files for duplication removal

    Optional parameters:
    --output_dir    Output directory for the processed BAM files and metrics (default: ./04_Markduplicates)
    --help          Show this help message
    """
    exit 0
}

// Validate required parameters
if (!params.input_dir) {
    exit 1, "ERROR: --input_dir must be specified"
}

// Pipeline default parameters
params.output_dir = params.get('output_dir', "${params.output_root}/04_Markduplicates")  // Output directory for processed files

// Mark and remove duplicated reads from the sorted BAM file 
process PICARD_MARKDUP {
    tag "Markduplicates on .sort.BAM: ${bam_file.simpleName}"

    input: 
        path bam_file
    
    output: 
        path "${params.output_dir}/${bam_file.simpleName}.dedup.bam"
        path "${params.output_dir}/${bam_file.simpleName}.metrics.txt"

    script:
    """
    picard MarkDuplicates \\
        --READ_NAME_REGEX null \\
        --INPUT ${bam_file} \\
        --OUTPUT ${params.output_dir}/${bam_file.simpleName}.dedup.bam \\
        --REMOVE_DUPLICATES \\
        --REMOVE_SEQUENCING_DUPLICATES \\
        --METRICS_FILE ${params.output_dir}/${bam_file.simpleName}.metrics.txt
    """
}

workflow {
    // Retrieve sorted BAM files from input directory
    Channel.fromPath("${params.input_dir}/*.sort.bam")
        | PICARD_MARKDUP()
        | view //debug
}