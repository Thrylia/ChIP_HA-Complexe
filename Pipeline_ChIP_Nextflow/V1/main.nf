#!/usr/bin/env nextflow

include { reportFastqc } from "./modules/fastqc.nf"
include { reportMultiqc } from "./modules/multiqc.nf"
include { cutadaptTrim } from "./modules/cutadapt.nf"

/*
    Default parameters
*/
params.batch = "Test"
params.adapter_file = "./data/adapters.fasta"
params.threads = "2"
params.length = "75"
params.quality = "20"

workflow {
    // FastQC report, WORKS
    samples_ch = Channel.fromPath("./data/samples.tsv")
                        .splitCsv(sep: '\t')
                        .map { r1, r2 -> [file(r1), file(r2)] } 
        //DEBUG : samples_ch.view { pair -> "R1: ${pair[0]} | R2: ${pair[1]}" }
    reportFastqc(samples_ch, params.batch)

    // MultiQC report
    reportMultiqc(reportFastqc.out.folderOut, params.batch)

    // Cutadapt trimming 
    cutadaptTrim(samples_ch, params.adapter_file, params.batch, params.threads, params.length, params.quality)

}
