#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN INFORMATION LOGGING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def runMessage() {
    log.info """
                              GENECE WISHBONE PIPELINE INFORMATION
        ===============================================================================

        INPUT OPTIONS
        -------------
        Input TSV:          ${params.input}
        Profiles:           ${workflow.profile}

        OUTPUT OPTIONS
        --------------
        Output Directory:   ${params.outdir}

        GC CORRECTION OPTIONS
        ---------------------
        GC Correction?      ${params.gc_corrected}

        GCPARAGON OPTIONS
        -----------------
        Genome:             ${params.rgb}


        SKIP OPTIONS
        ------------
        Skip FEMS?          ${params.skip_fems}
        Skip Coverage?      ${params.skip_coverage}


        PIPELINE VERSION:   ${params.version}

        ===============================================================================
        """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WISHBONE } from './workflows/wishbone'

//
// WORKFLOW: Run main genece/tfbscov pipeline
//
workflow GENECE_WISHBONE {
    runMessage()

    WISHBONE ()
}

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    GENECE_WISHBONE ()
}
