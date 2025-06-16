#!/usr/bin/env nextflow

/*
 * Trim the PE fastq files listed in an input_fastq file 
 * according to parameters given as arguments and a file 
 * containing the Illumina adapters used during sequencing,
 * and create a MultiQC report for the new reads.
 */

nextflow.enable.dsl=2

// Display help message
if (params.help) {
    log.info """
    Usage: nextflow run cutadapt_fastqc.nf --fastq_list <path> --adapter_file <path> [--length <int>] [--quality <int>] [--output_cutadapt <path>] [--output_trimqc <path>] [--help]

    Required parameters:
    --fastq_list    Path to a file containing the path to the paired-end FASTQ files to be trimmed
    --adapter_file  Path to a file containing the Illumina adapter sequences (Fasta file)

    Optional parameters:
    --length          Minimum length for the reads after trimming (default: 75)
    --quality         Minimum quality score for the reads after trimming (default: 20)
    --output_cutadapt Output directory for the trimmed FASTQ files (default: ./02_Trimmed-reads)
    --output_trimqc   Output directory for the FastQC reports (default: ./02_Trimmed-reads/FastQC_report)
    --help            Show this help message
    """
    exit 0
}

// Validate required parameters
if (!params.adapter_file) {
    exit 1, "ERROR: --adapter_file must be specified"
}
if (!params.fastq_list) {
    exit 1, "ERROR: --fastq_list must be specified"
}

// Pipeline default parameters
params.length = params.get('length', 75)                                                                // Length minimum for the reads
params.quality = params.get('quality', 20)                                                              // Quality minimum for the reads
params.output_cutadapt = params.get('output_cutadapt', "${params.output_root}/02_Trimmed-reads")                // Cutadapt output directory
params.output_trimqc = params.get('output_trimqc', "${params.output_root}/02_Trimmed-reads/FastQC_report")      // MultiQC output directory

// Cutadapt trimming
process TRIM_READS {
    tag "Trimming paired-end FASTQ files: ${replicate_id}"
    
    input: 
        tuple val(replicate_id), path(reads1), path(reads2)         // Tuple with replicate_id, Illumina R1 and R2
        path adapter_file = params.adapter_file                     // Illumina adapters

    output:
        tuple val(replicate_id),
        path("${params.output_cutadapt}/${replicate_id}_R1_trimmed.fastq.gz"),
        path("${params.output_cutadapt}/${replicate_id}_R2_trimmed.fastq.gz") into trimmed_reads
    
    script: 
    """
    cutadapt \\
        -a file:${adapter_file} -A file:${adapter_file} \\
        -o ${params.output_cutadapt}/${replicate_id}_R1_trimmed.fastq.gz \\
        -p ${params.output_cutadapt}/${replicate_id}_R2_trimmed.fastq.gz \\
        -j ${params.threads} \\
        --minimum-length ${params.length}:${params.length} \\
        -q ${params.quality} \\
        --pair-filter=any \\
        --quality-base=33 \\
        ${reads1} ${reads2}
    """
 }

// FastQC reports, quality check of the trimmed fastq files
process TRIM_QUALITY_CHECK {
    tag "FastQC of trimmed paired-end FASTQ files: ${replicate_id}"

    input: 
        tuple val(replicate_id), path(reads1_trimmed), path(reads2_trimmed) from trimmed_reads

    output:
        path "${params.output_trimqc}/data" into fastqc_reports
    
    script: 
    """
    mkdir -p "${params.output_trimqc}/data"
    fastqc -o "${params.output_trimqc}/data" -t ${params.threads} ${reads1_trimmed} ${reads2_trimmed}
    """
 }

// MultiQC report, it merges FastQC reports into one single report  
process TRIM_GROUP_QC {
    tag "Combining trimmed FastQC reports"
    
    input:
        path fastqc_reports_dir from fastqc_reports           // FastQC reports
    
    output:
        path "${params.output_trimqc}/multiqc_data"
        path "${params.output_trimqc}/multiqc_report.html"

    script: 
    """
    multiqc ${fastqc_reports_dir} -o "${params.output_trimqc}"
    """
}

workflow {
    // Create a channel from the file list containing paths to fastq files
    trimmed_reads = Channel.fromPath(params.fastq_list)
    | map { line -> 
        def fields = line.split('\t') 
        def replicate_id = fields[0].replaceAll(/_R1.fastq.gz$/, '') 
        tuple(replicate_id, file(fields[0]), file(fields[1])) 
    }
    | TRIM_READS()                                          // Stores the output in trimmed_reads

    fastqc_reports = trimmed_reads | TRIM_QUALITY_CHECK()   // Trimmed_reads as input + stores the output in fastqc_reports
    fastqc_reports | TRIM_GROUP_QC()                        // Fastqc_reports as input
 }