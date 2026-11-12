/*
 * Training workload for the AOT cache of the bundle (see bundle/train.sh).
 * It covers the common code paths of a run: config, includes, params, file
 * staging, operators, publishing and the local executor.
 */

include { COUNT } from './modules/count'

params.input = "${projectDir}/samplesheet.csv"
params.outdir = 'results'

process SPLIT {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.part'), emit: parts

    script:
    """
    awk 'NR % 2 == 1 { f = sprintf("${meta.id}_%d.part", ++n) } { print > f }' ${reads}
    """
}

workflow {
    main:
    samples = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> tuple([id: row.sample, group: row.group], file("${projectDir}/${row.reads}", checkIfExists: true)) }

    SPLIT(samples)
    COUNT(SPLIT.out.parts.transpose())

    counts = COUNT.out
        .map { meta, n -> tuple(meta.group, n as int) }
        .groupTuple()
        .map { group, ns -> "${group}\t${ns.sum()}" }
        .collectFile(name: 'counts.tsv', newLine: true, sort: true)

    counts.view { f -> f.text.trim() }

    publish:
    counts = counts
}

output {
    counts {
        path '.'
    }
}
