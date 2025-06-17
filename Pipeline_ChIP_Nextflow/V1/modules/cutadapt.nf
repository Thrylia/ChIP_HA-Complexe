#!/usr/bin/env nextflow


/*
    Run Cutadapt on the illumina.fastq.gz
    conda create -n cutadaptenv -c bioconda cutadapt=2.6
*/

process cutadaptTrim {
    publishDir 'results', mode: 'copy'
    
    input:
        tuple val (name), path (r1), path (r2)
        path adapter_file
        val batch
        val threads
        val length
        val quality

    output:
        path "${batch}_cutadaptTrim/${name}", emit : cutadaptDir
        path "${batch}_cutadaptTrim/${name}/TRIM_*.fastq.gz"

    script:
    """
    mkdir -p ${batch}_cutadaptTrim/${name}
    conda run -n cutadaptenv cutadapt \\
        -a file:${adapter_file} -A file:${adapter_file} \\
        -o ${batch}_cutadaptTrim/${name}/TRIM_${r1} \\
        -p ${batch}_cutadaptTrim/${name}/TRIM_${r2} \\
        -j ${threads} \\
        --minimum-length ${length}:${length} \\
        -q ${quality} \\
        --pair-filter=any \\
        --quality-base=33 \\
        ${r1} ${r2}
    """
}