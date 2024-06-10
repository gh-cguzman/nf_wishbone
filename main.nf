#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN INFORMATION LOGGING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// ANSII Colours used for terminal logging
//
public static Map logColours(Boolean monochrome_logs) {
    Map colorcodes = [:]

    // Reset / Meta
    colorcodes['reset']      = monochrome_logs ? '' : "\033[0m"
    colorcodes['bold']       = monochrome_logs ? '' : "\033[1m"
    colorcodes['dim']        = monochrome_logs ? '' : "\033[2m"
    colorcodes['underlined'] = monochrome_logs ? '' : "\033[4m"
    colorcodes['blink']      = monochrome_logs ? '' : "\033[5m"
    colorcodes['reverse']    = monochrome_logs ? '' : "\033[7m"
    colorcodes['hidden']     = monochrome_logs ? '' : "\033[8m"

    // Regular Colors
    colorcodes['black']      = monochrome_logs ? '' : "\033[0;30m"
    colorcodes['red']        = monochrome_logs ? '' : "\033[0;31m"
    colorcodes['green']      = monochrome_logs ? '' : "\033[0;32m"
    colorcodes['yellow']     = monochrome_logs ? '' : "\033[0;33m"
    colorcodes['blue']       = monochrome_logs ? '' : "\033[0;34m"
    colorcodes['purple']     = monochrome_logs ? '' : "\033[0;35m"
    colorcodes['cyan']       = monochrome_logs ? '' : "\033[0;36m"
    colorcodes['white']      = monochrome_logs ? '' : "\033[0;37m"

    // Bold
    colorcodes['bblack']     = monochrome_logs ? '' : "\033[1;30m"
    colorcodes['bred']       = monochrome_logs ? '' : "\033[1;31m"
    colorcodes['bgreen']     = monochrome_logs ? '' : "\033[1;32m"
    colorcodes['byellow']    = monochrome_logs ? '' : "\033[1;33m"
    colorcodes['bblue']      = monochrome_logs ? '' : "\033[1;34m"
    colorcodes['bpurple']    = monochrome_logs ? '' : "\033[1;35m"
    colorcodes['bcyan']      = monochrome_logs ? '' : "\033[1;36m"
    colorcodes['bwhite']     = monochrome_logs ? '' : "\033[1;37m"

    // Underline
    colorcodes['ublack']     = monochrome_logs ? '' : "\033[4;30m"
    colorcodes['ured']       = monochrome_logs ? '' : "\033[4;31m"
    colorcodes['ugreen']     = monochrome_logs ? '' : "\033[4;32m"
    colorcodes['uyellow']    = monochrome_logs ? '' : "\033[4;33m"
    colorcodes['ublue']      = monochrome_logs ? '' : "\033[4;34m"
    colorcodes['upurple']    = monochrome_logs ? '' : "\033[4;35m"
    colorcodes['ucyan']      = monochrome_logs ? '' : "\033[4;36m"
    colorcodes['uwhite']     = monochrome_logs ? '' : "\033[4;37m"

    // High Intensity
    colorcodes['iblack']     = monochrome_logs ? '' : "\033[0;90m"
    colorcodes['ired']       = monochrome_logs ? '' : "\033[0;91m"
    colorcodes['igreen']     = monochrome_logs ? '' : "\033[0;92m"
    colorcodes['iyellow']    = monochrome_logs ? '' : "\033[0;93m"
    colorcodes['iblue']      = monochrome_logs ? '' : "\033[0;94m"
    colorcodes['ipurple']    = monochrome_logs ? '' : "\033[0;95m"
    colorcodes['icyan']      = monochrome_logs ? '' : "\033[0;96m"
    colorcodes['iwhite']     = monochrome_logs ? '' : "\033[0;97m"

    // Bold High Intensity
    colorcodes['biblack']    = monochrome_logs ? '' : "\033[1;90m"
    colorcodes['bired']      = monochrome_logs ? '' : "\033[1;91m"
    colorcodes['bigreen']    = monochrome_logs ? '' : "\033[1;92m"
    colorcodes['biyellow']   = monochrome_logs ? '' : "\033[1;93m"
    colorcodes['biblue']     = monochrome_logs ? '' : "\033[1;94m"
    colorcodes['bipurple']   = monochrome_logs ? '' : "\033[1;95m"
    colorcodes['bicyan']     = monochrome_logs ? '' : "\033[1;96m"
    colorcodes['biwhite']    = monochrome_logs ? '' : "\033[1;97m"

    return colorcodes
}

public static String dashedLine(monochrome_logs) {
    Map colors = logColours(monochrome_logs)
    return "-${colors.dim}----------------------------------------------------${colors.reset}-"
}

def runMessage() {

    def monochrome_logs = params.monochrome_logs
    
    Map colors = logColours(monochrome_logs)

    log.info """
${colors.bipurple}GENECE WISHBONE PIPELINE INFORMATION${colors.reset}
${dashedLine(monochrome_logs)}

${colors.green}INPUT OPTIONS${colors.reset}
${dashedLine(monochrome_logs)}
Input TSV:          ${params.input}
Input BAM:          ${params.bams}
Profiles:           ${workflow.profile}
${dashedLine(monochrome_logs)}

${colors.green}OUTPUT OPTIONS${colors.reset}
${dashedLine(monochrome_logs)}
Output Directory:   ${params.outdir}
${dashedLine(monochrome_logs)}

${colors.green}GC CORRECTION OPTIONS${colors.reset}
${dashedLine(monochrome_logs)}
GC Correction?      ${params.gc_correction}
${dashedLine(monochrome_logs)}

${colors.green}END MOTIF CORRECTION OPTIONS${colors.reset}
${dashedLine(monochrome_logs)}
EM Correction?      ${params.em_correction}
2bit File:          ${params.genome_2bit}
Normalization:      ${params.em_norm}
${dashedLine(monochrome_logs)}

${colors.green}GCPARAGON OPTIONS${colors.reset}
${dashedLine(monochrome_logs)}
Genome:             ${params.rgb}
${dashedLine(monochrome_logs)}

${colors.green}TFBSCOV OPTIONS${colors.reset}
${dashedLine(monochrome_logs)}
Save Matrix?        ${params.save_tfbscov_matrix}
Window:             ${params.window}
Batch Size:         ${params.batch_size}
${dashedLine(monochrome_logs)}

${colors.green}SKIP OPTIONS${colors.reset}
${dashedLine(monochrome_logs)}
Skip FEMS?          ${params.skip_fems}
Skip Coverage?      ${params.skip_coverage}
Skip TFBSCov?       ${params.skip_tfbscov}
${dashedLine(monochrome_logs)}

${colors.green}PIPELINE VERSION:${colors.reset}   ${params.version}
${dashedLine(monochrome_logs)}
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
