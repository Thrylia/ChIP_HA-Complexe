#!/usr/bin/env nextflow

/*
 * Downsamples the different replicates so that 
 * they all contain the same number of read pairs and 
 * merges them into a single large sorted BAM file.
 */

nextflow.enable.dsl=2

// Display help message
if (params.help) {
    log.info """
    Usage: nextflow run downsample_merge.nf --input_dir <path> [--output_down <path>] [--output_merge <path>] [--help]

    Required parameters:
    --input_dir      Directory containing the BAM files to process

    Optional parameters:
    --output_down    Output directory for downsampled BAM files (default: ./06_Downsampled)
    --output_merge   Output directory for merged BAM files (default: ./07_Merged)
    --help           Show this help message
    """
    exit 0
}

// Validate required parameters
if (!params.input_dir) {
    exit 1, "ERROR: --input_dir must be specified"
}

// Pipeline default parameters
params.output_down = params.get('output_down', "${params.output_root}/06_Downsampled")       // Output directory for downsampled BAM files
params.output_merge = params.get('output_merge', "${params.output_root}/07_Merged")           // Output directory for merged BAM files

// Groups BAM files by prefix
process GROUP_BAM {
    tag "Grouping BAM files by prefix"

    input:
    path bam_files from bam_files_channel

    output:
    tuple val(prefix), path(bam_files) into grouped_bam_channel

    script:
    // Defines group names according to the prefix preceding 'RepX.sort.bam'
    def groups = bam_files.groupBy { it.name.replaceFirst(/-Rep\d+\.sort\.bam$/, '') }
    // Groups files by group name 
    def grouped_files = groups.collect { prefix, files ->
        tuple(prefix, files)
    }
    return grouped_files
}

// Counts the read pairs in BAM files within a group and returns the minimum found
process COUNT_READ_PAIRS {
    tag "Counting read pairs in BAM files"

    input:
    tuple val(prefix), path(bam_files)

    output:
    tuple val(prefix), path("${prefix}_counts.txt"), path(bam_files)

    script:
    """
    for bam in ${bam_files.join(' ')}; do
        samtools view -c -f 1 $bam > ${bam}.count
    done
    cat ${bam_files.collect { it + ".count" }.join(' ')} | sort -n | head -1 > ${prefix}_counts.txt
    """
}

// Downsamples the BAM files so that they all have the same number of read pairs
process SUBSAMPLE_BAM {
    tag "Subsampling BAM files"

    input:
    tuple val(prefix), path(count_file), path(bam_files)

    output:
    tuple val(prefix), path(bam_files.collect { "${params.output_down}/" + it.baseName + "_subsampled.bam" })

    script:
    // Retrieves the minimum number of read pairs in the group
    def min_reads = new File(count_file).text.trim()
    """
    for bam in ${bam_files.join(' ')}; do
        samtools view -s 42.${min_reads} -f 1 -b $bam > ${params.output_down}/${bam.baseName}_subsampled.bam
    done
    """
}

// Merges BAM files together, ensuring read pairs are kept together
process MERGE_BAM {
    tag "Merging BAM files"

    input:
    tuple val(prefix), path(bam_files)

    output:
    path "${params.output_merge}/${prefix}.merged.bam"

    script:
    def output_bam = "${params.output_merge}/${prefix}.merged.bam"
    def input_bams = bam_files.collect { it.toString() }.join(' ')
    """
    samtools merge -f ${output_bam} ${input_bams}
    """
}

// Sorts the merged BAM files by position
process SORT_BAM {
    tag "Sorting merged BAM files"

    input:
    path bam_file

    output:
    path "${params.output_merge}/${bam_file.baseName.replace('.bam', '.sorted.bam')}"

    script:
    """
    samtools sort -o ${params.output_merge}/${bam_file.baseName.replace('.bam', '.sorted.bam')} $bam_file
    
    # Remove the original (non-sorted) BAM file after sorting
    rm -f $bam_file
    """
}

workflow {
    // Create a channel from the BAM files in the input directory
    bam_files_channel = Channel.fromPath("${params.input_dir}/*.bam")
    
    // Process the BAM files: grouping, counting read pairs, downsampling, merging, and sorting
    grouped_bam_channel = GROUP_BAM(bam_files_channel)
    counted_bam_channel = COUNT_READ_PAIRS(grouped_bam_channel)
    subsampled_bam_channel = SUBSAMPLE_BAM(counted_bam_channel)
    merged_bam_channel = MERGE_BAM(subsampled_bam_channel)
    sorted_bam_channel = SORT_BAM(merged_bam_channel)
}