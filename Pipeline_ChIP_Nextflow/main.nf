#!/usr/bin/env nextflow

include { ALIGN_AND_FILTER } from "modules/bowtie_samtools.nf"
include { TRIM_READS; TRIM_QUALITY_CHECK; TRIM_GROUP_QC } from "modules/cutadapt_fastqc.nf"
include { COVERAGE_NORM; MATRIX_COMPUTE; GENERATE_PLOT } from "modules/deeptools.nf"
include { GROUP_BAM; COUNT_READ_PAIRS; SUBSAMPLE_BAM; MERGE_BAM; SORT_BAM } from "modules/downsample_merge.nf"
include { QUALITY_CHECK; GROUP_QC } from "modules/fastqc.nf"
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
//Deeptools, do I include it ?
if (!params.check_regions) {
    exit 1, "ERROR: --check_regions must be specified"
}
if (!(params.ref_point in ['center', 'TSS'])) {
    exit 1, "ERROR: --ref_point must be 'center' or 'TSS'. Provided: ${params.ref_point}"
}
// Outputs
params.output_fastqc = params.get('output_fastqc', "${params.output_root}/01_FastQC-report/data")           // FastQC Output directory
params.output_multiqc = params.get('out_multiqc', "${params.output_root}/01_FastQC-report")      // MultiQC Output directory                                                          // Quality minimum for the reads
params.output_cutadapt = params.get('output_cutadapt', "${params.output_root}/02_Trimmed-reads")                          // Cutadapt Output directory
params.output_trimqc = params.get('output_trimqc', "${params.output_root}/02_Trimmed-reads/FastQC_report")      // MultiQC Output directory
params.output_bowtie = params.get('output_bowtie', "${params.output_root}/03_Alignment")           // Bowtie2 Output directory
params.output_picard = params.get('output_picard', "${params.output_root}/04_Markduplicates")  // Output directory for processed files
params.output_intersect = params.get('output_intersect', "${params.output_root}/05_Blacklisted-Regions-Removed")    // Output directory
params.output_deeptools = params.get('output_deeptools', "${params.output_deeptools}/10_Deeptools")
params.output_down = params.get('output_down', "${params.output_root}/06_Downsampled")       // Output directory for downsampled BAM files
params.output_merge = params.get('output_merge', "${params.output_root}/07_Merged")           // Output directory for merged BAM files
params.output_peak = params.get('output_peak', "${params.output_root}/08_MACS-peaks")
params.output_annot = params.get('output_annot', "${params.output_root}/09_Annotated-peaks")
// Parameters for analyses
// Cutadapt
params.length = params.get('length', 75)                                                                // Length minimum for the reads
params.quality = params.get('quality', 20)    
// Deeptools
params.effective_size = params.get('effective_size', 2495461690)
params.length_min = params.get('length_min', 80)
params.length_max = params.get('length_max', 200)
params.smooth = params.get('smooth', 80)
params.before = params.get('before', 1000)
params.after = params.get('after', 1000)
params.plot_labels = params.get('plot_labels', '')
params.plot_title = params.get('plot_title', '')
// MACS
params.warning_file = params.get('warning_file', "${params.output_peak}/warnings.txt")
params.pe_mode = params.get('pe_mode', true)

 workflow {
    
    //get the variable value from the command line, create the option --greeting
    greeting_ch = Channel.of(params.greeting)
    
    // launch the process defined above
    name(greeting_ch)

 }