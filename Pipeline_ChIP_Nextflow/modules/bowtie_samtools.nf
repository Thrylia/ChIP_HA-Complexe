#!/usr/bin/env nextflow

/*
 * Sensitively aligns with Bowtie2 the reads stored 
 * in the input folder after creating an index of 
 * the reference genome, and outputs a sorted BAM file.
 */

nextflow.enable.dsl=2

// Display help message
if (params.help) {
    log.info """
    Usage: nextflow run bowtie.nf --input_dir <path> --genome <path> [--output_dir <path>]

    Required parameters:
    --input_dir    Path to the input directory containing FASTQ files
    --genome       Path to the reference genome file

    Optional parameters:
    --output_dir   Path to the output directory (default: ./03_Alignment)
    --help         Show this help message
    """
    exit 0
}

// Validate required parameters
if (!params.input_dir) {
    exit 1, "ERROR: --input_dir must be specified"
}
if (!params.genome) {
    exit 1, "ERROR: --genome must be specified"
}

// Pipeline default parameters
params.output_dir = params.get('output_dir', "${params.output_root}/03_Alignment")           // Bowtie2 Output directory

// Retrieve fastq files and group paired-end reads
Channel
    .fromPath("${params.input_dir}/*.fastq.gz")
    .map { file -> 
        def match = file.name =~ /(.+)_R[12]_trimmed.fastq.gz/
        if (match) {
            return tuple(file, match[0][1])     // Extract sample ID
        }
        return null
    }
    .filter { it != null }                      // Remove null entries
    .groupTuple()
    .set { fastq_pairs }

// Bowtie2 Alignment process
process ALIGN_AND_FILTER {
    tag "Aligning trimmed reads: ${replicate_id}"

    input:
        tuple path(reads1), path(reads2), val(replicate_id) from fastq_pairs

    output:
        path("${params.output_dir}/${replicate_id}.sort.bam"),
        path("${params.output_dir}/${replicate_id}.raw.bam"),
        path("${params.genome}.index")   // The index path is output for the next run, if required

    script:
    """
    # Build Bowtie2 index if not already present
    if [ ! -f "${params.genome}.index.1.bt2" ]; then
        bowtie2-build \\
            ${params.genome} ${params.genome}.index \\
            --threads ${params.threads}
    fi

    # Perform Bowtie2 alignment
    bowtie2 -x ${params.genome}.index \\
        -1 ${reads1} -2 ${reads2} \\
        -S ${replicate_id}.sam \\
        --threads ${params.threads} \\
        --phred33 \\
        --no-mixed \\
        --no-discordant \\
        --dovetail \\
        --very-sensitive

    # Convert SAM to BAM
    samtools view \\
        -bS ${replicate_id}.sam > ${replicate_id}.raw.bam

    # Convert SAM to sorted BAM + filter alignments
    samtools view \\
        -f 2 \\
        -F 1804 \\
        -q 20 \\
        -bS ${replicate_id}.sam > ${replicate_id}.tmp.bam
    samtools sort -O BAM \\
        -o ${replicate_id}.sort.bam ${replicate_id}.tmp.bam
    
    # Cleanup intermediate files
    rm ${replicate_id}.sam ${replicate_id}.tmp.bam
    """
}

// --very-sensitive => very-sensitive and not only sensitive
// --phred33 => to be sure that phred score taken is 33
// --no-mixed => only consider alignment status of pairs per se
// --no-discordant => the readings are directed towards each other
// --dovetail => the reads can overlap themselves
// --no-unal => don't output in SAM the unaligned PE, REMOVED FOR STATISTICS PURPOSE

// -q 20 : have mapping quality >= 20 , https://doi-org.insb.bib.cnrs.fr/10.1038/s41467-022-35384-1
// -F 1804 : Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates 
// -F 1284 : Remove non-primary alignments, duplicates and unmapped reads
// -f 2 : Retain properly paired reads

workflow {
    ALIGN_AND_FILTER()
}