nextflow.enable.types = true

include { MINI } from './workflows/mini'

params {
    // Samplesheet: sample,fastq_1,fastq_2[,clip_r1][,<tool>_args...]
    input: Path

    fasta: Path

    rrbs: Boolean

    // Default clip for R1; overridable per sample via `clip_r1` column
    clip_r1: Integer = 0

    clip_r2: Integer = 0

    skip_trimming_presets: Boolean

    skip_dedup: Boolean

    maxins: Integer?

    multiqc_title: String?

    // Per-tool args overrides, e.g. --args.trimgalore '--quality 30'
    args: Map<String,String> = [:]
}

workflow {
    main:
    ch_samples = channel.fromList(params.input.splitCsv(header: true)).map { row ->
        record(
            id: row.sample,
            single_end: !row.fastq_2,
            reads: row.fastq_2 ? [file(row.fastq_1), file(row.fastq_2)] : [file(row.fastq_1)],
            clip_r1: row.clip_r1 ? row.clip_r1.toInteger() : params.clip_r1,
            tool_args: row
                .findAll { e -> e.key.endsWith('_args') && e.value }
                .collectEntries { e -> [e.key.replace('_args', ''), e.value] }
        )
    }

    mini = MINI(ch_samples, params.fasta, params)

    publish:
    samples = mini.results
    multiqc = mini.multiqc_report
}

output {
    samples {
        path { r ->
            r.trim_log >> 'trimgalore/'
            r.align_log >> 'align/'
            r.sorted_bam >> 'align/'
            r.dedup_log >> 'dedup/'
            r.dedup_sorted_bam >> 'dedup/'
        }
        index {
            path 'samples.json'
        }
    }

    multiqc {
        path 'multiqc'
    }
}
