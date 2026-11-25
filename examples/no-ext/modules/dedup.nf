nextflow.enable.types = true

process DEDUP {
    tag id

    input:
    record(
        id: String,
        bam: Path,
        args: String?,
        prefix: String?
    )

    output:
    record(
        id: id,
        dedup_bam: file('*.dedup.bam'),
        dedup_log: file('*.log')
    )

    script:
    def pfx = prefix ?: id
    """
    echo "dedup ${args ?: ''} ${bam}" > ${pfx}.dedup.log
    touch ${pfx}.dedup.bam
    """
}
