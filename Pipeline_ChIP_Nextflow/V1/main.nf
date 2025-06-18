#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//nextflow run main.nf -resume -profile docker,test

include { reportRawFastqc; reportTrimFastqc } from "./modules/fastqc.nf"
include { reportRawMultiqc; reportTrimMultiqc } from "./modules/multiqc.nf"
include { cutadaptTrim } from "./modules/cutadapt.nf"
include { bowtieIndex; bowtieAlign } from "./modules/bowtie.nf"
include { samtoolsFilter; samtoolsSubset; samtoolsSort; samtoolsSortWhite } from "./modules/samtools.nf"
include { picardValidate; picardReplace; picardDuplicates } from "./modules/picard.nf"
include { bedtoolsBlackListed } from "./modules/bedtools.nf"

workflow {

    main:

    // Raw samples, tuple(name, r1, r2)
    samples_ch = Channel.fromPath(params.samples)
                        .splitCsv(sep: '\t')
                        //.view(pair -> "name: ${pair[0]} | R1: ${pair[1]} | R2: ${pair[2]}")

    // FastQC report
    reportRawFastqc(samples_ch, params.batch)

    // MultiQC report 
    reportRawMultiqc(reportRawFastqc.out.fastqcDir.collect(), params.batch)

    // Cutadapt trimming 
    cutadaptTrim(samples_ch, params.adapter_fasta, params.batch, params.threads, params.length, params.quality)
    def trimReads_ch = cutadaptTrim.out.cutadaptTup

    // FastQC report 
    reportTrimFastqc(trimReads_ch, params.batch)

    // MultiQC report 
    reportTrimMultiqc(reportTrimFastqc.out.fastqcDir.collect(), params.batch)

    // Bowtie2 index
    bowtieIndex(params.genome_fasta, params.batch)

    // Bowtie Align
    bowtieAlign(trimReads_ch,bowtieIndex.out.bowtieIndexDir, params.threads, params.batch)

    // Samtools filters + SAM to BAM, remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates
    samtoolsFilter(bowtieAlign.out.bowtieAligned, params.batch)

    // Picard Validate
    picardValidate(samtoolsFilter.out.bamSorted, params.genome_fasta, params.batch) 
    // If not validated
    // presence of a small pirouette, forcing picard to be true. To be modified in V2
    picardReplace(picardValidate.out.validateFile, samtoolsFilter.out.bamSorted, params.batch) 

    // Picard no Duplicates
    picardDuplicates(picardReplace.out.validateBam, params.batch)
    samtoolsSort(picardDuplicates.out.noDupBam, params.batch)

    // Bedtools, no black listed regions
    bedtoolsBlackListed(samtoolsSort.out.noDupSortBam, params.blacklisted_regions, params.batch)
    samtoolsSortWhite(bedtoolsBlackListed.out.bamWhite, params.batch)
}
