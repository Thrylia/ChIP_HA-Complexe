#!/usr/bin/env nextflow


/*
    Picard processes to MarkDuplicates
        git clone https://github.com/broadinstitute/picard.git
        cd picard/
        ./gradlew shadowJar
        java -jar ~/Softwares/picard/build/libs/picard.jar => launch picard 
*/


process picardValidate {
    //errorStrategy 'ignore'

    publishDir 'results', mode: 'copy'

    input:
        tuple val (name), path (sortBam)
        path genomeFasta
        val batch

    output:
        tuple val (name), path ("${batch}_6-checkBam/${name}/${name}.sort.bam.validate"), emit: validateFile

    script:
    """
    mkdir -p ${batch}_6-checkBam/${name}
    java -jar ~/Softwares/picard/build/libs/picard.jar ValidateSamFile \\
        -I ${sortBam} \\
        -M SUMMARY \\
        -O ${batch}_6-checkBam/${name}/${name}.sort.bam.validate \\
        -R ${genomeFasta} || true
    """
}

process picardReplace {
    publishDir 'results', mode: 'copy'
    
    input:
        tuple val (name), path (validateFile)
        tuple val (nameBam), path (sortBam)
        val batch

    output:
        tuple val (nameBam), path ("${batch}_5-filteredAlign/${nameBam}/${nameBam}.sort.bam"), emit: validateBam
        path "${batch}_5-filteredAlign/${nameBam}/${nameBam}.old.sort.bam", optional: true

    script:
    """
    mkdir -p ${batch}_5-filteredAlign/${nameBam}
    if [[ \$(cat ${validateFile} | grep "ERROR:MISSING_READ_GROUP") ]]; then
        echo "Validation report exists, proceeding with AddOrReplaceReadGroups"
        java -jar ~/Softwares/picard/build/libs/picard.jar AddOrReplaceReadGroups \\
            -I ${sortBam} \\
            -O ${batch}_5-filteredAlign/${nameBam}/${nameBam}.tmp.bam \\
            -RGLB lib2024 \\
            -RGPL illumina \\
            -RGPU run1 \\
            -RGSM ${nameBam} || true

        mv ${sortBam} ${batch}_5-filteredAlign/${nameBam}/${nameBam}.old.sort.bam 
        mv ${batch}_5-filteredAlign/${nameBam}/${nameBam}.tmp.bam ${batch}_5-filteredAlign/${nameBam}/${nameBam}.sort.bam
    else
        echo "Skipping ${nameBam}: no valid .validate file or empty"
        # Output for Nextflow :
        cp $sortBam ${batch}_5-filteredAlign/${nameBam}/${nameBam}.sort.bam
    fi
    """
}

process picardDuplicates {
    publishDir 'results', mode: 'copy'

    input:
        tuple val (name), path (validateBam)
        val batch

    output:
        tuple val (name), path ("${batch}_7-noDuplicates/${name}/${name}.sort.bam"), emit: noDupBam
        path "${batch}_7-noDuplicates/${name}/${name}.picstats"

    script:
    """
    mkdir -p ${batch}_7-noDuplicates/${name}
    java -jar ~/Softwares/picard/build/libs/picard.jar MarkDuplicates \\
        --INPUT ${validateBam} \\
        --OUTPUT ${batch}_7-noDuplicates/${name}/${name}.bam \\
        --REMOVE_DUPLICATES \\
        --METRICS_FILE ${batch}_7-noDuplicates/${name}/${name}.picstats 
    samtools sort -O BAM \\
        -o ${batch}_7-noDuplicates/${name}/${name}.sort.bam \\
        ${batch}_7-noDuplicates/${name}/${name}.bam
    rm ${batch}_7-noDuplicates/${name}/${name}.bam
    """
}
