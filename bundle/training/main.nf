/*
 * Training workload for the AOT cache of the bundle (see bundle/train.sh).
 * It covers the common code paths of a run with static typing: config,
 * includes, typed params and samplesheet loading, records, file staging,
 * operators, exec tasks, publishing and the local executor.
 */
nextflow.enable.types = true

include { COUNT } from './modules/count'

params {
    input: Channel<Sample>
    outdir: String = 'results'
}

record Sample {
    id: String
    group: String
    reads: Path
}

process SPLIT {
    tag id

    input:
    record(
        id: String,
        group: String,
        reads: Path
    )

    output:
    record(
        id: id,
        group: group,
        parts: files('*.part')
    )

    script:
    """
    awk 'NR % 2 == 1 { f = sprintf("${id}_%d.part", ++n) } { print > f }' ${reads}
    """
}

process SUMMARY {
    input:
    name: String
    lines: List<String>

    output:
    file(name)

    exec:
    def path = task.workDir.resolve(name)
    path.text = lines.join('\n') + '\n'
}

workflow {
    main:
    parts = SPLIT(params.input)
        .flatMap { s -> s.parts.collect { part -> record(id: s.id, group: s.group, part: part) } }

    lines = COUNT(parts)
        .map { r -> tuple(r.group, r.count.trim().toInteger()) }
        .groupBy()
        .map { group, ns -> "${group}\t${ns.sum()}" }
        .collect()
        .map { ls -> ls.toSorted() }

    counts = SUMMARY('counts.tsv', lines)
    counts.view { f -> f.text.trim() }

    publish:
    counts = counts
}

output {
    counts: Path {
        path '.'
    }
}
