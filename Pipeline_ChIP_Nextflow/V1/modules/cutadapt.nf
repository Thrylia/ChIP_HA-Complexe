#!/usr/bin/env nextflow


/*
    Run Cutadapt on the illumina.fastq.gz
    conda create -n cutadaptenv -c bioconda cutadapt=2.6
*/

process cutadaptTrim {
    publishDir 'results', mode: 'copy'
    
    input:
        tuple path (R1), path (R2)
        path adapter_file
        val batch
        val threads
        val length
        val quality

    output:
        path "${batch}_cutadaptTrim", emit : folderOut
        path "${batch}_cutadaptTrim/TRIM_*.fastq.gz"

    script:
    """
    mkdir ${batch}_cutadaptTrim
    conda run -n cutadaptenv cutadapt \\
        -a file:${adapter_file} -A file:${adapter_file} \\
        -o ${batch}_cutadaptTrim/TRIM_${R1} \\
        -p ${batch}_cutadaptTrim/TRIM_${R2} \\
        -j ${threads} \\
        --minimum-length ${length}:${length} \\
        -q ${quality} \\
        --pair-filter=any \\
        --quality-base=33 \\
        ${R1} ${R2}
    """
}