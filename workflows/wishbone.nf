/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ONE       } from '../modules/local/samtools/index/main.nf'
include { GCPARAGON                                  } from '../modules/local/gcparagon/main.nf'
include { CREATE_FEMS_MATRIX                         } from '../modules/local/custom/fems/main.nf'
include { CREATE_COVERAGE_MATRIX                     } from '../modules/local/custom/coverage/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow WISHBONE {

    runMessage()

    ch_versions = Channel.empty()

    // Create a channel from input
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

    //
    // MODULE: CREATE INDEX FOR EACH BAM USING SAMTOOLS
    //
    SAMTOOLS_INDEX_ONE( ch_input )

    //
    // MODULE: CORRECT FOR GC BIAS USING GCPARAGON
    //
    if (!params.gc_corrected) {
        GCPARAGON( SAMTOOLS_INDEX_ONE.out.bam_bai )

        ch_bam_bai = GCPARAGON.out.bam_bai
    } else {
        ch_bam_bai = SAMTOOLS_INDEX_ONE.out.bam_bai
    }

    //
    // MODULE: CREATE FEMS MATRIX
    //
    if (!params.skip_fems) {
        CREATE_FEMS_MATRIX( ch_bam_bai )
    }

    //
    // MODULE: CREATE COVERAGE MATRIX
    //
    if (!params.skip_coverage) {
        CREATE_COVERAGE_MATRIX(
            ch_bam_bai,
            file(params.regions),
            file(params.blacklist)
        )
    }
}