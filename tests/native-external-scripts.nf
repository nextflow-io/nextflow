#!/usr/bin/env nextflow

/*
 * Native external scripts contract fixture.
 *
 * This file intentionally uses the proposed ADR syntax from
 * adr/20260522-native-external-scripts.md. It is listed in
 * tests/checks/.IGNORE until native external scripts are implemented.
 */

params.reads = "$baseDir/data/sample.fa"
params.prefix = 'native'

process bash_sort {
    cpus 2

    input:
    path reads

    output:
    path 'bash.sorted.fa', emit: sorted
    path 'bash.context.txt', emit: context

    script:
    file 'scripts/native-external/bash_sort.sh'
}

process python_summary {
    input:
    path reads
    val prefix

    output:
    path 'python.summary.txt', emit: summary
    path 'python.context.json', emit: context

    script:
    file 'scripts/native-external/python_summary.py'
}

process r_summary {
    input:
    path reads
    val prefix

    output:
    path 'r.summary.txt', emit: summary
    path 'r.context.txt', emit: context

    script:
    file 'scripts/native-external/r_summary.R'
}

workflow {
    reads_ch = Channel.of(file(params.reads))

    bash_sort(reads_ch)
    python_summary(reads_ch, params.prefix)
    r_summary(reads_ch, params.prefix)

    bash_sort.out.context.view { it.text }
    python_summary.out.summary.view { it.text }
    r_summary.out.summary.view { it.text }
}
