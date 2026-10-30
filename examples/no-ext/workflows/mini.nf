nextflow.enable.types = true

include { TRIMGALORE                          } from '../modules/trimgalore'
include { ALIGN                               } from '../modules/align'
include { DEDUP                               } from '../modules/dedup'
include { SAMTOOLS_SORT                       } from '../modules/samtools_sort'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_DEDUP } from '../modules/samtools_sort'
include { MULTIQC                             } from '../modules/multiqc'
include { Sample ; MiniParams                 } from './types'
include { toolArgs ; trimgaloreArgs ; alignArgs } from './args'

workflow MINI {
    take:
    ch_samples: Channel<Sample>
    fasta: Path
    params: MiniParams

    main:
    ch_trimmed = TRIMGALORE(
        ch_samples.map { s ->
            s + record(args: toolArgs('trimgalore', s, params, trimgaloreArgs(s, params)))
        }
    )

    // re-join by id to recover sample fields (single_end etc.)
    ch_aligned = ALIGN(
        ch_samples.join(ch_trimmed, by: 'id').map { r ->
            r + record(reads: r.trim_reads, fasta: fasta, args: toolArgs('align', r, params, alignArgs(r, params)))
        }
    )

    ch_sorted = SAMTOOLS_SORT(
        ch_samples.join(ch_aligned, by: 'id').map { r ->
            record(id: r.id, bam: r.bam, args: toolArgs('samtools_sort', r, params, ''), prefix: "${r.id}.sorted")
        }
    )

    def ch_dedup: Channel<DedupResult> = channel.empty()
    if( !params.skip_dedup ) {
        ch_dedup_raw = DEDUP(
            ch_samples.join(ch_aligned, by: 'id').map { r ->
                record(id: r.id, bam: r.bam, args: toolArgs('dedup', r, params, ''))
            }
        )
        // same module, different call site -> different prefix (was withName selector)
        ch_dedup_sorted = SAMTOOLS_SORT_DEDUP(
            ch_dedup_raw.map { r -> record(id: r.id, bam: r.dedup_bam, prefix: "${r.id}.dedup.sorted") }
        )
        ch_dedup = ch_dedup_raw
            .join(ch_dedup_sorted.map { r -> record(id: r.id, dedup_sorted_bam: r.sorted_bam) }, by: 'id')
            .map { r -> record(id: r.id, dedup_log: r.dedup_log, dedup_sorted_bam: r.dedup_sorted_bam) }
    }

    ch_results = ch_samples
        .join(ch_trimmed, by: 'id')
        .join(ch_aligned, by: 'id')
        .join(ch_sorted, by: 'id')
        .join(ch_dedup, by: 'id', remainder: true)
        .map { r ->
            record(
                id: r.id,
                single_end: r.single_end,
                trim_log: r.trim_log,
                align_log: r.align_log,
                sorted_bam: r.sorted_bam,
                dedup_log: r.dedup_log,
                dedup_sorted_bam: r.dedup_sorted_bam
            )
        }

    ch_logs = ch_results
        .flatMap { r -> [r.trim_log, r.align_log] }
        .collect()
    multiqc_report = MULTIQC(ch_logs, params.args['multiqc'] ?: (params.multiqc_title ? "--title \"${params.multiqc_title}\"" : null))

    emit:
    results: Channel<SampleResult> = ch_results
    multiqc_report: Value<Path> = multiqc_report
}

record DedupResult {
    id: String
    dedup_log: Path
    dedup_sorted_bam: Path
}

record SampleResult {
    id: String
    single_end: Boolean
    trim_log: Path
    align_log: Path
    sorted_bam: Path
    dedup_log: Path?
    dedup_sorted_bam: Path?
}
