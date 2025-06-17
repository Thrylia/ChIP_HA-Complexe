#!/usr/bin/env nextflow


/*
    Run Cutadapt on the illumina.fastq.gz
    conda create -n cutadaptenv -c bioconda cutadapt=2.6
*/

process cutadaptTrim {
    
    input:
        tuple val (name), path (r1), path (r2)
        path adapter_file
        val batch
        val threads
        val length
        val quality

    output:
    tuple val(name), path("${batch}_2-cutadaptTrim/${name}/TRIM_${r1}"), path("${batch}_2-cutadaptTrim/${name}/TRIM_${r2}"), emit: cutadaptTup
    path "${batch}_2-cutadaptTrim/${name}/TRIM_*.fastq.gz"

    script:
    """
    mkdir -p ${batch}_2-cutadaptTrim/${name}
    conda run -n cutadaptenv cutadapt \\
        -a file:${adapter_file} -A file:${adapter_file} \\
        -o ${batch}_2-cutadaptTrim/${name}/TRIM_${r1} \\
        -p ${batch}_2-cutadaptTrim/${name}/TRIM_${r2} \\
        -j ${threads} \\
        --minimum-length ${length}:${length} \\
        -q ${quality} \\
        --pair-filter=any \\
        --quality-base=33 \\
        ${r1} ${r2}
    """
}