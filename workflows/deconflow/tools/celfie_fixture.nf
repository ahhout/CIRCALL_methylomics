#!/usr/bin/env nextflow
// Exercise the pinned native process; the script is supplied by upstream bin/.
include { RUN_DECONV } from "${params.deconflow_source}/modules/celfie/deconv/main.nf"

workflow {
    RUN_DECONV(file(params.count_table), 2)
}
