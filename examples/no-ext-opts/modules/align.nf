nextflow.enable.types = true

process ALIGN {
    tag id

    input:
    record(
        id: String,
        single_end: Boolean,
        reads: List<Path>,
        fasta: Path,
        args: String?,
        prefix: String?
    )

    output:
    record(
        id: id,
        bam: file('*.bam'),
        align_log: file('*.log')
    )

    script:
    def pfx = prefix ?: id
    def fastq = single_end ? reads.join(' ') : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    echo "align ${args ?: ''} --genome ${fasta} ${fastq}" > ${pfx}.align.log
    touch ${pfx}.bam
    """
}
