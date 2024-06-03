/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ONE       } from '../modules/local/samtools/index/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TWO       } from '../modules/local/samtools/index/main.nf'
include { GCPARAGON                                  } from '../modules/local/gcparagon/main.nf'
include { EMCORRECTION                               } from '../modules/local/emcorrection/main.nf'
include { CREATE_FEMS_MATRIX                         } from '../modules/local/custom/fems/main.nf'
include { CREATE_COVERAGE_MATRIX                     } from '../modules/local/custom/coverage/main.nf'
include { CREATE_TFBSCOV_MATRIX                      } from '../modules/local/custom/tfbscov/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow WISHBONE {

    ch_versions = Channel.empty()

    ch_2bit = Channel.fromPath(params.genome_2bit, checkIfExists: true)
    ch_blacklist = Channel.fromPath(params.blacklist, checkIfExists: true)
    ch_regions = Channel.fromPath(params.regions, checkIfExists: true)

    // Create a channel from input
    if (params.input) {
        Channel
            .fromPath( params.input )
            .splitCsv( header: true, sep: '\t' )
            .map {
                row ->
                    def meta = [:]
                    meta.sample_id = row.sample_id
                    [ meta, file(row.bam_path) ]
            }
            .set { ch_input }
    } else if (params.bams) {
        Channel
            .fromPath(params.bams, checkIfExists: true)
            .map {
                bam ->
                    def meta = [:]
                    meta.sample_id = bam.baseName
                    [ meta, bam ]
            }
            .set { ch_input }
    } else if (params.input && params.bam) {
        error "You cannot set both --input and --bams! Choose one or the other."
    } else if (!params.input && !params.bam) {
        error "You must set one of --input or --bams! You set neither."
    }

    //
    // MODULE: CREATE INDEX FOR EACH BAM USING SAMTOOLS
    //
    SAMTOOLS_INDEX_ONE( ch_input )

    SAMTOOLS_INDEX_ONE
        .map {
            meta, bam, bai ->
                [ meta, bam, bai ]
        }
        .set { ch_samtools }

    //
    // MODULE: CORRECT FOR GC BIAS USING GCPARAGON
    //
    if (params.gc_correction) {
        GCPARAGON( ch_samtools )

        GCPARAGON
        .map {
            meta, bam, bai ->
                [ meta, bam, bai ]
        }
        .set { ch_gcparagon }

    } else {
        ch_samtools
        .map {
            meta, bam, bai ->
                [ meta, bam, bai ]
        }
        .set { ch_gcparagon }
    }

    //
    // MODULE: CORRECT FOR END MOTIF BIAS USING CUSTOM
    //
    if (params.em_correction) {
        EMCORRECTION(
            ch_gcparagon,
            ch_2bit,
            ch_blacklist
        )

        SAMTOOLS_INDEX_TWO( EMCORRECTION.out.bam )

        SAMTOOLS_INDEX_TWO
        .map {
            meta, bam, bai ->
                [ meta, bam, bai ]
        }
        .set { ch_emcorrection }

    } else {
        ch_gcparagon
        .map {
            meta, bam, bai ->
                [ meta, bam, bai ]
        }
        .set { ch_emcorrection }
    }

    //
    // MODULE: CREATE FEMS MATRIX
    //
    if (!params.skip_fems) {
        CREATE_FEMS_MATRIX( ch_emcorrection )
    }

    //
    // MODULE: CREATE COVERAGE MATRIX
    //
    if (!params.skip_coverage) {
        CREATE_COVERAGE_MATRIX(
            ch_emcorrection,
            ch_regions,
            ch_blacklist
        )
    }

    //
    // MODULE: CREATE TFBS FEATURES
    //
    if (!params.skip_tfbscov) {
        CREATE_TFBSCOV_MATRIX( ch_emcorrection )
    }
}