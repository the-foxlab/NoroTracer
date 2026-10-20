#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    the-foxlab/norotracer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/the-foxlab/norotracer
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NOROTRACER  } from './workflows/norotracer'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_norotracer_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_norotracer_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow THEFOXLAB_NOROTRACER {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    NOROTRACER (
        samplesheet,
        params.outdir,
    )

    emit:
    versions       = NOROTRACER.out.versions                 // channel: [ path(versions.yml) ]
    similarities_csv = NOROTRACER.out.similarities_csv
    reduced_genbank_file = NOROTRACER.out.reduced_genbank_file
    split_orf1_output = NOROTRACER.out.split_orf1_output
    split_orf2_output = NOROTRACER.out.split_orf2_output
    orf1_aln_output = NOROTRACER.out.orf1_aln_output
    orf2_aln_output = NOROTRACER.out.orf2_aln_output
    tree_output = NOROTRACER.out.tree_output
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    THEFOXLAB_NOROTRACER (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.monochrome_logs,
    )

    publish:
    similarities_csv = THEFOXLAB_NOROTRACER.out.similarities_csv
    reduced_genbank_file = THEFOXLAB_NOROTRACER.out.reduced_genbank_file
    split_orf1_output = THEFOXLAB_NOROTRACER.out.split_orf1_output
    orf1_aln_output = THEFOXLAB_NOROTRACER.out.orf1_aln_output
    orf2_aln_output = THEFOXLAB_NOROTRACER.out.orf2_aln_output
    tree_output = THEFOXLAB_NOROTRACER.out.tree_output
    split_orf2_output = THEFOXLAB_NOROTRACER.out.split_orf2_output
}

output {
    similarities_csv {
        path {meta, csv -> "similarity_outputs/${meta.id}_similarities/"} 
        mode 'copy'
    }
    reduced_genbank_file {
        path 'filtered_references'
        mode 'copy'
    }
    split_orf1_output {
        path {meta, fasta -> "split_orfs/${meta.id}_orf1/"}
        mode 'copy'
    }
    split_orf2_output {
        path {meta, fasta -> "split_orfs/${meta.id}_orf2/"}
        mode 'copy'
    }
    orf1_aln_output {
        path {aln -> "alignments/orf1/"}
        mode 'copy'
    }
    orf2_aln_output {
        path {aln -> "alignments/orf2/"}
        mode 'copy'
    }
    tree_output {
        path {meta, tree -> "trees/${meta.id}/"}
        mode 'copy'
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
