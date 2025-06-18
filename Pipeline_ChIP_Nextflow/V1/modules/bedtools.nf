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
        tuple val (name), path ("${batch}_8-noBlackListedRegions/${name}/${name}.bam"), emit: bamWhite

    script:
    """
    mkdir -p ${batch}_8-noBlackListedRegions/${name}
    bedtools intersect \\
        -abam ${bamSorted} \\
        -b ${blackListedRegions} \\
        -v \\
        -sorted \\
        > ${batch}_8-noBlackListedRegions/${name}/${name}.bam
    """
}