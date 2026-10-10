nextflow.enable.types = true

process TRIMGALORE {
    tag id

    input:
    record(
        id: String,
        single_end: Boolean,
        reads: List<Path>,
        args: String?,
        prefix: String?
    )

    output:
    record(
        id: id,
        trim_reads: files('*_trimmed.fq.gz').toSorted(),
        trim_log: file('*.log')
    )

    script:
    def pfx = prefix ?: id
    def trimmed = single_end ? "${pfx}_trimmed.fq.gz" : "${pfx}_1_trimmed.fq.gz ${pfx}_2_trimmed.fq.gz"
    """
    echo "trim_galore ${args ?: ''} ${single_end ? '' : '--paired'} ${reads.join(' ')}" > ${pfx}.log
    touch ${trimmed}
    """
}
