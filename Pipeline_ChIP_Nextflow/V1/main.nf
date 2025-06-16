#!/usr/bin/env nextflow

include { reportFastqc } from "/home/apeticca/Documents/Test_NextFlow/modules/fastqc.nf"
include { reportMultiqc } from "/home/apeticca/Documents/Test_NextFlow/modules/multiqc.nf"

/*
    Default parameters
*/
params.batch = "Test"


workflow {
    // FastQC report, WORKS
    samples_ch = Channel.fromPath("/home/apeticca/Documents/Test_NextFlow/data/samples.tsv")
                        .splitCsv(sep: '\t')
                        .map { r1, r2 -> [file(r1), file(r2)] } 
        //DEBUG : samples_ch.view { pair -> "R1: ${pair[0]} | R2: ${pair[1]}" }
    reportFastqc(samples_ch, params.batch)

    // MultiQC report
    reportMultiqc(reportFastqc.out.folderOut, params.batch)

}