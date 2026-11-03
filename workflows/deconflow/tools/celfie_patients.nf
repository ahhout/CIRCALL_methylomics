#!/usr/bin/env nextflow
// Invoke the unchanged pinned native fitting process with explicitly prepared counts.
include { RUN_DECONV } from "${params.deconflow_source}/modules/celfie/deconv/main.nf"

process OBSERVE_RESTARTS {
    label 'process_low'
    input:
    path count_table
    path diagnostic_script
    path native_script
    val nsamples
    output:
    path 'restart_diagnostics.json', emit: diagnostics
    path '*tissue_proportions.txt', emit: fractions
    script:
    """
    python ${diagnostic_script} ${native_script} restart_diagnostics.json ${count_table} diagnostic ${nsamples} -m ${params.celfie_maxiter} -u ${params.unknown} -c ${params.converg} -r ${params.celfie_randrest} -p ${params.parall_job}
    """
}

workflow {
    RUN_DECONV(file(params.count_table), params.nsamples)
    if (params.diagnostics) {
        OBSERVE_RESTARTS(file(params.count_table), file(params.diagnostic_script),
                         file("${params.deconflow_source}/bin/celfie.py"), params.nsamples)
    }
}
