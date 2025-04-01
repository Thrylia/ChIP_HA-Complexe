#!/usr/bin/env nextflow

include { QUALITY_CHECK; GROUP_QC } from "modules/fastqc.nf"
include { ALIGN_AND_FILTER } from "modules/bowtie_samtools.nf"
include { TRIM_READS; TRIM_QUALITY_CHECK; TRIM_GROUP_QC } from "modules/cutadapt_fastqc.nf"
include { COVERAGE_NORM; MATRIX_COMPUTE; GENERATE_PLOT } from "modules/deeptools.nf"
include { GROUP_BAM; COUNT_READ_PAIRS; SUBSAMPLE_BAM; MERGE_BAM; SORT_BAM } from "modules/downsample_merge.nf"
include { REMOVE_REGIONS; GET_STATS } from "modules/intersect_stats.nf"
include { GET_BAM_FILES; PEAKS_CALLING; ANNOTATE_PEAKS } from "modules/macs_homer.nf"
include { PICARD_MARKDUP } from "modules/picard_markduplicates.nf"


// Mandatory paramaters 
if (!params.fastq_list) {
    exit 1, "ERROR: --fastq_list must be specified"
}
if (!params.adapter_file) {
    exit 1, "ERROR: --adapter_file must be specified"
}
if (!params.input_dir) {
    exit 1, "ERROR: --input_dir must be specified"
}
if (!params.genome) {
    exit 1, "ERROR: --genome must be specified"
}
if (!params.gtf_file) {
    error "ERROR: --gtf_file must be specified"
}
if (!params.black_regions) {
    exit 1, "ERROR: --black_regions must be specified"
}
// !!!!!! Deeptools, not included in the main. Too many regions to check, best to do it with the module directly !!!!!!
/*
if (!params.check_regions) {
    exit 1, "ERROR: --check_regions must be specified"
}
if (!(params.ref_point in ['center', 'TSS'])) {
    exit 1, "ERROR: --ref_point must be 'center' or 'TSS'. Provided: ${params.ref_point}"
}
*/
// !!!!!! Deeptools, not included in the main. Too many regions to check, best to do it with the module directly !!!!!!

// Outputs
params.output_fastqc = params.get('output_fastqc', "${params.output_root}/01_FastQC-report/data")                   // FastQC Output directory
params.output_multiqc = params.get('out_multiqc', "${params.output_root}/01_FastQC-report")                         // MultiQC Output directory                                                          // Quality minimum for the reads
params.output_cutadapt = params.get('output_cutadapt', "${params.output_root}/02_Trimmed-reads")                    // Cutadapt Output directory
params.output_trimqc = params.get('output_trimqc', "${params.output_root}/02_Trimmed-reads/FastQC_report")          // MultiQC Output directory
params.output_bowtie = params.get('output_bowtie', "${params.output_root}/03_Alignment")                            // Bowtie2 Output directory
params.output_picard = params.get('output_picard', "${params.output_root}/04_Markduplicates")                       // Markduplicates Output directory 
params.output_intersect = params.get('output_intersect', "${params.output_root}/05_Blacklisted-Regions-Removed")    // bedtools intersect Output directory  
params.output_down = params.get('output_down', "${params.output_root}/06_Downsampled")                              // downsampled BAM files Output directory 
params.output_merge = params.get('output_merge', "${params.output_root}/07_Merged")                                 // merged BAM files Output directory
params.output_peak = params.get('output_peak', "${params.output_root}/08_MACS-peaks")                               // MACS Output directory 
params.output_annot = params.get('output_annot', "${params.output_root}/09_Annotated-peaks")                        // Homer Output directory 
params.output_deeptools = params.get('output_deeptools', "${params.output_deeptools}/10_Deeptools")                 // deepTools Output directory 

// Parameters for analyses
// Cutadapt
params.length = params.get('length', 75)                                                    // Length filter reads
params.quality = params.get('quality', 20)                                                  // Quality filter reads
// MACS
params.warning_file = params.get('warning_file', "${params.output_peak}/warnings.txt")      // Warning file name
params.pe_mode = params.get('pe_mode', true)                                                // Data paired-end

// !!!!!! Deeptools, not included in the main. Too many regions to check, best to do it with the module directly !!!!!!
// Deeptools
/*
params.effective_size = params.get('effective_size', 2495461690)                            // Genome effective size
params.length_min = params.get('length_min', 80)                                            // Minimum length reads
params.length_max = params.get('length_max', 200)                                           // Maximum length reads
params.smooth = params.get('smooth', 80)                                                    // Plot smoothing
params.before = params.get('before', 1000)                                                  // Bases before regions to plot
params.after = params.get('after', 1000)                                                    // Bases after regions to plot
params.plot_labels = params.get('plot_labels', '')                                          // Samples labels, plot legend
params.plot_title = params.get('plot_title', '')                                            // Plot title
*/
// !!!!!! Deeptools, not included in the main. Too many regions to check, best to do it with the module directly !!!!!!

 workflow {
    // Quality report of the raw reads
    fastqc_reports = QUALITY_CHECK(params.fastq_list)
    GROUP_QC(fastqc_reports)
    // Trimming
    trimmed_reads = Channel.fromPath(params.fastq_list)
    | map { line -> 
        def fields = line.split('\t') 
        def replicate_id = fields[0].replaceAll(/_R1.fastq.gz$/, '') 
        tuple(replicate_id, file(fields[0]), file(fields[1])) 
    }
    TRIM_READS(trimmed_reads)
    // Quality report of the trimmed reads
    fastqc_reports = TRIM_QUALITY_CHECK(trimmed_reads)
    TRIM_GROUP_QC(fastqc_reports)
    // Align to genome
    sorted_bam = ALIGN_AND_FILTER(trimmed_reads).out.filter { it.endsWith('.sort.bam') }
    // Remove the PCR duplicates
    dedup_bam = PICARD_MARKDUP(sorted_bam).out.filter { it.endsWith('.dedup.bam') }
    // Remove the black listed regions
    white_bam = REMOVE_REGIONS(dedup_bam)
    GET_STATS(white_bam)
    // Downsamples
    sorted_bam = GROUP_BAM(white_bam) | COUNT_READ_PAIRS() | SUBSAMPLE_BAM () | MERGE_BAM () | SORT_BAM()
    // MACS
    query_control_pairs = GET_BAM_FILES(params.output_merge)
    query_control_pairs.view().each { pair ->
        def bam_file_query = pair[0]
        def bam_file_control = pair[1]
        def prefix_out = bam_file_query.name.replace(".sort.bam", "")
    narrowPeak_file = PEAKS_CALLING(bam_file_query, bam_file_control, params.pe_mode, prefix_out).out.filter { it.name == "${prefix_out}_peaks.narrowPeak.bed" }.first()
    ANNOTATE_PEAKS(narrowPeak_file, prefix_out, params.genome_file, params.gtf_file)
 }