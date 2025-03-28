#!/usr/bin/env nextflow

/*
 * Removes ENCODE-blacklisted regions from the alignment and
 * sorts the new BAM file created, followed by generating stats.
 */

nextflow.enable.dsl=2

// Display help message
if (params.help) {
    log.info """
    Usage: nextflow run intersect_stats.nf --input_dir <path> --black_regions <path> [--output <path>] [--help]

    Required parameters:
    --input_dir      Directory containing the BAM files to process
    --black_regions  Path to a file containing the blacklisted regions (BED format)

    Optional parameters:
    --output         Output directory for the processed BAM files (default: ./05_Blacklisted-Regions-Removed)
    --help           Show this help message
    """
    exit 0
}

// Validate required parameters
if (!params.input_dir) {
    exit 1, "ERROR: --input_dir must be specified"
}
if (!params.black_regions) {
    exit 1, "ERROR: --black_regions must be specified"
}

// Pipeline default parameters
params.output_intersect_intersect = params.get('output_intersect', "${params.output_root}/05_Blacklisted-Regions-Removed")    // Output directory

// Remove the ENCODE-blacklisted regions from the alignment + Sort
process REMOVE_REGIONS {
    tag "Remove regions in *.dedup.bam: ${replicate_id}"

    input: 
        path bam_file
        path black_regions = params.black_regions

    output:
        path "${params.output_intersect}/${replicate_id}.sort.bam"
    
    script:
    def replicate_id = bam_file.baseName.replaceFirst(/\.dedup\.bam$/, '')  // Get replicate ID
    """
    # Remove the regions in the BAM file from the black_regions file 
    bedtools intersect \\
        -abam ${bam_file} \\
        -b ${black_regions} \\
        -v \\
        -sorted \\
        > ${params.output_intersect}/${replicate_id}.bam

    # Sort the new BAM file
    samtools sort \\
        -O BAM \\
        -o ${params.output_intersect}/${replicate_id}.sort.bam \\
        ${params.output_intersect}/${replicate_id}.bam

    # Cleanup the temporary file
    rm ${params.output_intersect}/${replicate_id}.bam
    """
}

// Generate statistics reports
process GET_STATS{
    tag "Generate stats for BAM: ${sorted_bam_file}"

    input:
        path sorted_bam_file

    output:
        path "${sorted_bam_file}.bai", "${sorted_bam_file}.flagstats", "${sorted_bam_file}.stats"

    script:
    """
    # Index the BAM file
    if [ ! -f ${sorted_bam_file}.bai ]; then
        samtools index -o ${sorted_bam_file}.bai -@ ${params.threads} ${sorted_bam_file}
    fi

    # Generate flagstat summary
    samtools flagstat -@ ${params.threads} -O tsv ${sorted_bam_file} > ${sorted_bam_file}.flagstats

    # Generate samtools stats report 
    samtools stats -@ ${params.threads} ${sorted_bam_file} > ${sorted_bam_file}.stats
    """
}

workflow {
    // Gets the BAM files from input directory
    list_bams = Channel.fromPath("${params.input_dir}/*.dedup.bam")

    // Process each BAM file
    processed_bams = list_bams | REMOVE_REGIONS()

    // Generate statistics for the processed BAM files
    processed_bams | GET_STATS()
}