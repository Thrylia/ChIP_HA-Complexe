#!/usr/bin/env nextflow

/*
    Remove blacklisted regions with bedtools
*/

process bedtoolsBlackListed {
    input:
        tuple val (name), path (bamSorted)
        path blackListedRegions
        val batch

    output:
        tuple val (name), path ("${batch}_8-noBlackListedRegions/${name}/${name}.sort.bam"), emit: bamWhite

    script:
    """
    mkdir -p ${batch}_8-noBlackListedRegions/${name}
    bedtools intersect \\
        -abam ${bamSorted} \\
        -b ${blackListedRegions} \\
        -v \\
        -sorted \\
        > ${batch}_8-noBlackListedRegions/${name}/${name}.bam
    samtools sort \\
        -O BAM \\
        -o ${batch}_8-noBlackListedRegions/${name}/${name}.sort.bam \\
        ${batch}_8-noBlackListedRegions/${name}/${name}.bam
    """
}



