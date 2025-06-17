#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { reportRawFastqc; reportTrimFastqc } from "./modules/fastqc.nf"
include { reportRawMultiqc; reportTrimMultiqc } from "./modules/multiqc.nf"
include { cutadaptTrim } from "./modules/cutadapt.nf"
include { bowtieIndex; bowtieAlign } from "./modules/bowtie.nf"

/*
    Default parameters
*/
params.batch = "Test"
params.adapter_fasta = "./data/adapters.fasta"
params.threads = "2"
params.length = "75"
params.quality = "20"
params.genome_fasta = "./data/Mus_musculus.GRCm39.chr1.fa"

workflow {
// Raw samples, tuple(name, r1, r2)
    samples_ch = Channel.fromPath("./data/samples.tsv")
                        .splitCsv(sep: '\t')
                        //.view(pair -> "name: ${pair[0]} | R1: ${pair[1]} | R2: ${pair[2]}")

    // FastQC report
    reportRawFastqc(samples_ch, params.batch)

    // MultiQC report 
    reportRawMultiqc(reportRawFastqc.out.fastqcDir.collect(), params.batch)

    // Cutadapt trimming 
    cutadaptTrim(samples_ch, params.adapter_fasta, params.batch, params.threads, params.length, params.quality)
    trimReads_ch = cutadaptTrim.out.cutadaptDir.collect().flatten()

    // FastQC report 
    reportTrimFastqc(trimReads_ch, params.batch)

    // MultiQC report 
    reportTrimMultiqc(reportTrimFastqc.out.fastqcDir.collect(), params.batch)

    // Bowtie2 index
    bowtieIndex(params.genome_fasta, params.batch)

    // Bowtie Align
    bowtieAlign(trimReads_ch,bowtieIndex.out.indexDir, params.threads, params.batch)
}
