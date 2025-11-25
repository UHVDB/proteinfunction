#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UHVDB/proteinfunction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/UHVDB/proteinfunction
----------------------------------------------------------------------------------------
    Overview:
        1. Download latest ICTV VMR (Nextflow)

*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
process HASHDEREP {
    label "process_low"

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("${meta.id}.fna.gz")                      , emit: faa
    tuple val(meta), path(".command.log")                           , emit: log
    tuple val(meta), path(".command.sh")                            , emit: script

    script:
    """
    # remove partial proteins
    seqkit grep \\
        ${faa} \\
        --pattern ".*partial=00.*" \\
        --by-name \\
        --use-regexp \\
        --threads ${task.cpus} \\
        --out-file ${meta.id}.full.faa

    # extract hash and length for each protein
    seqkit fx2tab \\
        ${meta.id}.full.faa \\
        --name \\
        --only-id \\
        --seq-hash \\
        --length \\
        --threads ${task.cpus} \\
        --out-file ${meta.id}.fx2tab.tsv

    # identify unique hashes + lengths
    csvtk uniq \\
        ${meta.id}.fx2tab.tsv \\
        --no-header-row \\
        --tabs \\
        --num-cpus ${task.cpus} \\
        --fields 2,3 | \\
    csvtk cut \\
        --no-header-row \\
        --tabs \\
        --fields 1 \\
        --num-cpus ${task.cpus} \\
        --out-file ${meta.id}.uniq.tsv
        
    # extract unique sequences
    seqkit grep \\
        ${meta.id}.full.faa \\
        --pattern-file ${meta.id}.uniq.tsv \\
        --threads ${task.cpus} \\
        --out-file ${meta.id}.uniq.faa.gz
    """
}


// Run entry workflow
workflow {

    main:
    // Check if output file already exists
    def output_file = file("${params.output}")
    def ch_faa = channel.fromPath(params.input_faa).map { faa ->
            [ [ id: "${faa.getBaseName()}" ], faa ]
        }

    if (!output_file.exists()) {

        // Remove partial proteins and dereplicate by hash+length
        HASHDEREP(
            ch_faa
        )


    } else {
        println "[UHVDB/proteinfunction]: Output file [${params.output}] already exists!"
    }

    // Delete intermediate and Nextflow-specific files
    def remove_tmp = params.remove_tmp
    workflow.onComplete {
        if (output_file.exists()) {
            def work_dir = new File("./work/")
            def nextflow_dir = new File("./.nextflow/")
            def launch_dir = new File(".")

            work_dir.deleteDir()
            nextflow_dir.deleteDir()
            launch_dir.eachFileRecurse { file ->
                if (file.name ==~ /\.nextflow\.log.*/) {
                    file.delete()
                }
            }

            if (remove_tmp) {
                def tmp_dir = new File("./tmp/")
                tmp_dir.deleteDir()
            }
        }
    }
}