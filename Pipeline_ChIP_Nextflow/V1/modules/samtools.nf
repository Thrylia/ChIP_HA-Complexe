#!/usr/bin/env nextflow


/*
    samtools processes
*/

process samtoolsFilter {
    input:
        tuple val (name), path (rawSam)
        val batch

    output:
        tuple val (name), path ("${batch}_5-filteredAlign/${name}/${name}.sort.bam"), emit: bamSorted

    script:
    """
    mkdir -p ${batch}_5-filteredAlign/${name}
    samtools view \\
        -F 1804 \\
        -f 2 \\
        -q 20 \\
        -bS ${rawSam} \\
        > ${batch}_5-filteredAlign/${name}/${name}.bam 
    samtools sort -O BAM \\
        -o ${batch}_5-filteredAlign/${name}/${name}.sort.bam \\
        ${batch}_5-filteredAlign/${name}/${name}.bam
    rm ${batch}_5-filteredAlign/${name}/${name}.bam
    """
}

process samtoolsSort {
    publishDir "results/${batch}_7-noDuplicates/${name}", mode: 'copy'
    
    input:
        tuple val (name), path (bamUnsorted)
        val batch

    output:
        tuple val (name), path ("${name}.sort.bam"), emit: noDupSortBam

    script:
    """
    samtools sort -O BAM \\
        -o ${name}.sort.bam \\
        ${bamUnsorted}
    rm ${bamUnsorted}
    """
}

process samtoolsSortWhite {    
    input:
        tuple val (name), path (bamUnsorted)
        val batch

    output:
        tuple val (name), path ("${batch}_8-noBlackListedRegions/${name}/${name}.sort.bam"), emit: whiteSortBam

    script:
    """
    mkdir -p ${batch}_8-noBlackListedRegions/${name}
    samtools sort -O BAM \\
        -o ${batch}_8-noBlackListedRegions/${name}/${name}.sort.bam \\
        ${bamUnsorted}
    rm ${bamUnsorted}
    """
}

process samtoolsSubset { 
    input:
        path filteredBamDir
        val batch

    output:
        path "\${batch}_Subset/*/*sort.bam"

    script:
    """
    fraction_to_keep=\$(samtools idxstats \$name.sort.bam \
        | cut -f3 | awk -v ct=\$number_reads_to_downsample 'BEGIN {total=0} {total += \$1} END {print ct/total}')
    samtools view \\
        -b \\
        --subsample \$fraction_to_keep \\
        -@ \$SLURM_CPUS_PER_TASK \\
        > \$sample_name.bam
    samtools sort \\
        -O BAM \\
        -o \$sample_name.sort.bam \\
        \$sample_name.bam
    """
}