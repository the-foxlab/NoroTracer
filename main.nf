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
    all_similarities_csv = NOROTRACER.out.all_similarities_csv
    reduced_genbank_file = NOROTRACER.out.reduced_genbank_file
    split_orf_output = NOROTRACER.out.split_orf_output
    comb_aln_output = NOROTRACER.out.comb_aln_output
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
    all_similarities_csv = THEFOXLAB_NOROTRACER.out.all_similarities_csv
    reduced_genbank_file = THEFOXLAB_NOROTRACER.out.reduced_genbank_file
    split_orf_output = THEFOXLAB_NOROTRACER.out.split_orf_output
    comb_aln_output = THEFOXLAB_NOROTRACER.out.comb_aln_output
    tree_output = THEFOXLAB_NOROTRACER.out.tree_output
}

output {
    similarities_csv {
        path {meta, csv -> 
            csv >> "similarity_outputs/per_sample/${meta.id}_similarities.${csv.extension}"
        } 
        mode 'copy'
    }
    all_similarities_csv {
        path {"similarity_outputs/"} 
        mode 'copy'
    }
    reduced_genbank_file {
        path 'filtered_references'
        mode 'copy'
    }
    split_orf_output {
        path {meta, fasta -> 
                fasta >> "split_orfs/${meta.orf}/${meta.id}_${fasta.baseName}.${fasta.extension}"            
        }
        mode 'copy'
    }
    comb_aln_output {
        path {meta, fasta -> "alignments/${meta.id}/"}
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
