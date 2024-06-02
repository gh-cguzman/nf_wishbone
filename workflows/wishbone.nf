/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ONE       } from '../modules/local/samtools/index/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TWO       } from '../modules/local/samtools/index/main.nf'
include { GCPARAGON                                  } from '../modules/local/gcparagon/main.nf'
include { CREATE_FEMS_MATRIX                         } from '../modules/local/custom/fems/main.nf'
include { CREATE_COVERAGE_MATRIX                     } from '../modules/local/custom/coverage/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow WISHBONE {

    ch_versions = Channel.empty()

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
        ch_input = Channel.fromPath(params.bams, checkIfExists: true)
    }

    //
    // MODULE: CREATE INDEX FOR EACH BAM USING SAMTOOLS
    //
    SAMTOOLS_INDEX_ONE( ch_input )

    //
    // MODULE: CORRECT FOR GC BIAS USING GCPARAGON
    //
    if (params.gc_correction) {
        GCPARAGON( SAMTOOLS_INDEX_ONE.out.bam_bai )

        ch_bam_bai_gc = GCPARAGON.out.bam_bai
    } else {
        ch_bam_bai_gc = SAMTOOLS_INDEX_ONE.out.bam_bai
    }

    //
    // MODULE: CORRECT FOR END MOTIF BIAS USING CUSTOM
    //
    if (params.gc_correction && params.em_correction) {
        EMCORRECTION(
            ch_bam_bai_gc,
            file(params.genome_2bit),
            file(params.blacklist)
        )

        SAMTOOLS_INDEX_TWO( EMCORRECTION.out.bam )

        ch_bam_bai_em = SAMTOOLS_INDEX_TWO.out.bam_bai
    } else {
        ch_bam_bai_em = ch_bam_bai_gc
    }

    //
    // MODULE: CREATE FEMS MATRIX
    //
    if (!params.skip_fems) {
        CREATE_FEMS_MATRIX( ch_bam_bai_em )
    }

    //
    // MODULE: CREATE COVERAGE MATRIX
    //
    if (!params.skip_coverage) {
        CREATE_COVERAGE_MATRIX(
            ch_bam_bai_em,
            file(params.regions),
            file(params.blacklist)
        )
    }
}