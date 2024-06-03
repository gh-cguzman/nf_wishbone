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

    ch_2bit = file(params.genome_2bit)
    ch_blacklist = file(params.blacklist)
    ch_regions = file(params.regions)

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

    //
    // MODULE: CORRECT FOR GC AND END MOTIF BIAS
    //
    if (params.gc_correction && params.em_correction) {
        GCPARAGON( SAMTOOLS_INDEX_ONE.out.bam_bai )

        EMCORRECTION(
            GCPARAGON.out.bam_bai,
            ch_2bit,
            ch_blacklist
        )

        SAMTOOLS_INDEX_TWO( EMCORRECTION.out.bam )

        ch_corrected_bams = SAMTOOLS_INDEX_TWO.out.bam_bai

    } else if (params.gc_correction && !params.em_correction) {
        GCPARAGON( SAMTOOLS_INDEX_ONE.out.bam_bai )

        ch_corrected_bams = GCPARAGON.out.bam_bai
    } else if (!params.gc_correction && params.em_correction) {
        EMCORRECTION(
            SAMTOOLS_INDEX_ONE.out.bam_bai,
            ch_2bit,
            ch_blacklist
        )

        SAMTOOLS_INDEX_TWO( EMCORRECTION.out.bam )

        ch_corrected_bams = SAMTOOLS_INDEX_TWO.out.bam_bai
    } else {
        error "Something wierd happened. Cry and contact Carlos."
    }

    //
    // MODULE: CREATE FEMS MATRIX
    //
    if (!params.skip_fems) {
        CREATE_FEMS_MATRIX( ch_corrected_bams )
    }

    //
    // MODULE: CREATE COVERAGE MATRIX
    //
    if (!params.skip_coverage) {
        CREATE_COVERAGE_MATRIX(
            ch_corrected_bams,
            ch_regions,
            ch_blacklist
        )
    }

    //
    // MODULE: CREATE TFBS FEATURES
    //
    if (!params.skip_tfbscov) {
        CREATE_TFBSCOV_MATRIX( ch_corrected_bams )
    }
}