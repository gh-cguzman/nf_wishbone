#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN INFORMATION LOGGING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def runMessage() {
    log.info """

                       RUN INFORMATION
        =================================================

        INPUT OPTIONS
        =============
        Input File:         ${params.input}
        Regions   :         ${params.regions}
        Blacklist :         ${params.blacklist}

        OUTPUT OPTIONS
        ==============
        Output Directory:   ${params.outdir}

        GC CORRECTION OPTIONS
        =====================
        GC Correction?      ${params.gc_corrected}

        SKIP OPTIONS
        ============
        Skip FEMS?          ${params.skip_fems}
        Skip Coverage?      ${params.skip_coverage}


        PIPELINE VERSION:    ${params.version}

        =================================================
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
    WISHBONE ()
}

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    GENECE_WISHBONE ()
}
