#!/usr/bin/env nextflow
// Input adaptation only: the native UXM fitting process is imported unchanged.
include { UXM_DECONV } from "${params.deconflow_source}/modules/uxm/uxm_deconv/main.nf"

// DecoNFlow does not emit counts. Use the same counter and l4 defaults as UXM.
process FRAGMENT_COUNTS {
    label 'process_low'
    input:
    tuple val(sample), path(pat), path(csi)
    path markers
    output:
    path "${sample}.uxm.bed.gz"
    script:
    """
    wgbstools homog --rlen ${params.fragment_rlen} -b ${markers} ${pat} --prefix ${sample}
    """
}

workflow {
    if (params.fragment_method == 'counts') {
        pats = Channel.fromPath(params.fragment_samples).splitCsv(header: true)
        FRAGMENT_COUNTS(pats.map { tuple(it.sample, file(it.pat), file(it.csi)) }, file(params.fragment_reference))
    } else if (params.fragment_method == 'uxm') {
        pats = Channel.fromPath(params.fragment_samples).splitCsv(header: true)
        UXM_DECONV(pats.map { file(it.pat) }.collect(),
                   pats.map { file(it.csi) }.collect(), file(params.fragment_reference))
    } else {
        error "Unsupported fragment method: ${params.fragment_method}; EpiDISH is retired from this entry point"
    }
}
