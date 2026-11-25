nextflow.enable.types = true

process SAMTOOLS_SORT {
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
        sorted_bam: file("${prefix ?: "${id}.sorted"}.bam")
    )

    script:
    """
    echo "samtools sort ${args ?: ''} ${bam}" > ${prefix ?: "${id}.sorted"}.bam
    """
}
