#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

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
