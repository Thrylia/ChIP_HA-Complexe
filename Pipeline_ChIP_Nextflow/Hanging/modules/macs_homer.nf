#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Display help message
if (params.help) {
    log.info """
    Usage: nextflow run macs_homer.nf --input_dir <path> --genome_file <path> --gtf_file <path> [--output_peak <path>] [--output_annot <path>] [--warning_file <path>] [--pe_mode <true/false>]

    Required parameters:
    --input_dir    Path to the input directory containing BAM files
    --genome_file  Path to the genome file
    --gtf_file     Path to the GTF file

    Optional parameters:
    --output_peak  Output directory for MACS peak files (default: ./08_MACS-peaks)
    --output_annot Output directory for HOMER annotated peaks (default: ./09_Annotated-peaks)
    --warning_file File to store warnings (default: ./08_MACS-peaks/warnings.txt)
    --pe_mode      Enable paired-end mode (default: true)
    --help         Show this help message
    """
    exit 0
}

// Validate required parameters
if (!params.input_dir) {
    error "ERROR: --input_dir must be specified"
}
if (!params.genome_file) {
    error "ERROR: --genome_file must be specified"
}
if (!params.gtf_file) {
    error "ERROR: --gtf_file must be specified"
}

// Default pipeline parameters
params.output_peak = params.get('output_peak', "${params.output_root}/08_MACS-peaks")
params.output_annot = params.get('output_annot', "${params.output_root}/09_Annotated-peaks")
params.warning_file = params.get('warning_file', "${params.output_peak}/warnings.txt")
params.pe_mode = params.get('pe_mode', true)

// Retrieves BAM files and matches query/control pairs
process GET_BAM_FILES {
    input:
    path input_dir

    output:
    tuple path(bam_file_query), path(bam_file_control)

    script:
    """
    def bam_files = file(input_dir).listFiles().findAll { it.name.endsWith('.sort.bam') }
    def query_control_pairs = []

    bam_files.each { bam_file_query ->
        if (!bam_file_query.name.startsWith("Input-")) {
            def bam_file_control_name = "Input-" + bam_file_query.name
            def bam_file_control = bam_files.find { it.name == bam_file_control_name }

            if (bam_file_control) {
                query_control_pairs << [bam_file_query, bam_file_control]
            } else {
                new File("${params.warning_file}").append("WARNING: No control file found for ${bam_file_query.name}\n")
                query_control_pairs << [bam_file_query, null]
            }
        }
    }

    query_control_pairs.each { pair -> 
        emit(pair[0], pair[1])  // Directly emit the tuple here
    }
    """
}

// MACS3 peak calling for each query/control pair
process PEAKS_CALLING {
    input:
    path bam_file_query
    path bam_file_control
    val pe_mode
    val prefix_out

    output:
    path("${params.output_peak}/${prefix_out}*") into all_peaks_files

    script:
    """
    macs3 callpeak \\
        -t ${bam_file_query} \\
        ${bam_file_control ? "-c " + bam_file_control : "--nolambda"} \\
        -f ${pe_mode ? 'BAMPE' : 'BAM'} \\
        -g 2309746861 \\
        -n ${prefix_out} \\
        -B \\
        -q 0.01 \\
        --nomodel \\
        --shift 0 \\
        --extsize 147 \\
        --call-summits \\
        --outdir ${params.output_peak}
    """
}

// Annotation of peaks using HOMER
process ANNOTATE_PEAKS {
    input:
    path peaks_file
    val prefix_out
    val genome_file
    val gtf_file

    output:
    path("${params.output_annot}/${prefix_out}.annot.narrowPeaks.bed")

    script:
    """
    annotatePeaks.pl \\
        ${peaks_file} \\
        ${genome_file} \\
        -gtf ${gtf_file} \\
        > ${params.output_annot}/${prefix_out}.annot.narrowPeaks.bed
    """
}

workflow {
    // Retrieves BAM file query/control pairs
    query_control_pairs = GET_BAM_FILES(input_dir: params.input_dir)

    query_control_pairs.view().each { pair ->
        def bam_file_query = pair[0]
        def bam_file_control = pair[1]
        def prefix_out = bam_file_query.name.replace(".sort.bam", "")

        // Perform peak calling
        peaks_outputs = PEAKS_CALLING(bam_file_query, bam_file_control, params.pe_mode, prefix_out)

        // Select the correct *.narrowPeak.bed file for annotation
        narrowPeak_file = peaks_outputs.filter { it.name == "${prefix_out}_peaks.narrowPeak.bed" }.first()

        // Annotate peaks using HOMER
        ANNOTATE_PEAKS(narrowPeak_file, prefix_out, params.genome_file, params.gtf_file)
    }
}