#!/usr/bin/env nextflow

/*
 * Generates normalized BigWig files from BAM files,
 * Computes a matrix over a region of interest,
 * Creates an average signal profile plot (SVG).
 */

nextflow.enable.dsl = 2

// Display help message
if (params.help) {
    log.info """
    Usage: nextflow run script_name.nf --input_dir <path> --check_regions <path> --ref_point <point> [options]

    Required parameters:
    --input_dir    Directory containing BAM files to process
    --check_regions    File containing regions of interest for matrix computation
    --ref_point    Reference point for matrix computation, either 'center' or 'TSS'

    Optional parameters:
    --output       Output directory where results will be stored (default: <output_root>/10_Deeptools)
    --effective_size  Effective genome size for normalization (default: 2495461690)
    --length_min   Minimum fragment length for BAM file normalization (default: 80)
    --length_max   Maximum fragment length for BAM file normalization (default: 200)
    --smooth       Smoothing window size for BAM file normalization (default: 80)
    --before       Distance (in base pairs) to extend before region start for matrix computation (default: 1000)
    --after        Distance (in base pairs) to extend after region start for matrix computation (default: 1000)
    --plot_labels  Labels for the plot (comma-separated) (optional)
    --plot_title   Title for the plot (optional)
    --help         Show this help message
    """
    exit 0
}

// Validate required parameters
if (!params.input_dir) {
    exit 1, "ERROR: --input_dir must be specified"
}
if (!params.check_regions) {
    exit 1, "ERROR: --check_regions must be specified"
}
if (!(params.ref_point in ['center', 'TSS'])) {
    exit 1, "ERROR: --ref_point must be 'center' or 'TSS'. Provided: ${params.ref_point}"
}

// Define parameters with default values
params.output = params.get('output', "${params.output_root}/10_Deeptools")
params.effective_size = params.get('effective_size', 2495461690)
params.length_min = params.get('length_min', 80)
params.length_max = params.get('length_max', 200)
params.smooth = params.get('smooth', 80)
params.before = params.get('before', 1000)
params.after = params.get('after', 1000)
params.plot_labels = params.get('plot_labels', '')
params.plot_title = params.get('plot_title', '')

// Get BAM files to process (excluding 'Input-' files)
bam_files = file("${params.input_dir}/*merged.sort.bam").filter { !it.name.startsWith('Input-') }

// Generates normalized BigWig file from BAM file
process COVERAGE_NORM {
    tag "Normalizing BAM files: $sample_bw"

    input:
        path bam_file from bam_files

    output:
        path "${params.output}/${sample_bw}" into bw_outputs

    script:
    """
    def sample_bw = bam_file.baseName.replace('.merged.sort.bam', '.rpgc.bw')
    bamCoverage -b ${bam_file} \\
        --outFileFormat bigwig \\
        --outFileName ${params.output}/${sample_bw} \\
        --normalizeUsing RPGC \\
        --effectiveGenomeSize ${params.effective_size} \\
        --MNase \\
        --numberOfProcessors ${params.threads} \\
        --minFragmentLength ${params.length_min} \\
        --maxFragmentLength ${params.length_max} \\
        --smoothLength ${params.smooth}
    """
}

// Computes the matrix over a region of interest
process MATRIX_COMPUTE {
    tag "Computing matrix: $bw_file"

    input:
        path bw_file from bw_outputs

    output:
        path "${params.output}/${sample_id}.${regions_name}.rpgc.gz" into matrix_outputs

    script:
    """
    def sample_id = bw_file.baseName.replace('.rpgc.bw', '')
    def regions_name = params.check_regions.baseName
    computeMatrix reference-point \\
        -S ${bw_file} \\
        -R ${params.check_regions} \\
        --beforeRegionStartLength ${params.before} \\
        --afterRegionStartLength ${params.after} \\
        --referencePoint ${params.ref_point} \\
        --sortRegions keep \\
        --numberOfProcessors ${params.threads} \\
        --smartLabels \\
        --metagene \\
        --outFileName ${params.output}/${sample_id}.${regions_name}.rpgc.gz
    """
}

// Creates the average signal profile plot (SVG)
process GENERATE_PLOT {
    tag "Generating SVG plot: ${sample_id}"

    input:
        path matrix_file from matrix_outputs

    output:
        path "${params.output}/${sample_id}.profile.svg"

    script:
    """
    def plotTitleOption = params.plot_title ? "--plotTitle '${params.plot_title}'" : ""
    def samplesLabelOption = params.plot_labels ? "--samplesLabel '${params.plot_labels}'" : ""
    def sample_id = matrix_file.baseName.replace('.rpgc.gz', '')
    
    plotProfile -m ${matrix_file} \\
        -o ${params.output}/${sample_id}.profile.svg \\
        --plotFileFormat svg \\
        --perGroup \\
        --plotType lines \\
        --numPlotsPerRow 1 \\
        --plotHeight 9 \\
        --plotWidth 7 \\
        --legendLocation best \\
        ${plotTitleOption} \\
        ${samplesLabelOption}
    """
}

workflow {
    COVERAGE_NORM() | MATRIX_COMPUTE() | GENERATE_PLOT()
}