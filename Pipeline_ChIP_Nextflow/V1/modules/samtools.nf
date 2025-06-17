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


process samtoolsSubset {
    publishDir 'results', mode 'copy'

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



