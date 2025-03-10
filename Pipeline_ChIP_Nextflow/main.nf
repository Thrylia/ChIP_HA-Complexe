#!/usr/bin/env nextflow

/*
 * Pipeline default parameters
 */
params.greeting = "Holà mundo!"


/*
 * Command to print a chr in output.txt
 */
 process name {

    //outputDir to use
    publishDir 'results', mode: 'copy'

    input: 
        val greeting

    output:
        path 'output.txt'
    
    script: 
    """
    echo "$greeting" > output.txt
    """

 }


 workflow {
    
    //get the variable value from the command line, create the option --greeting
    greeting_ch = Channel.of(params.greeting)
    
    // launch the process defined above
    name(greeting_ch)

 }