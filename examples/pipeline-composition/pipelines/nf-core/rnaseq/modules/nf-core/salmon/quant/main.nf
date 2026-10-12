nextflow.enable.types = true

process SALMON_QUANT {
    tag "${id}"

    input:
    record(
        id: String,
        bam: Path
    )

    output:
    record(
        id: id,
        counts: file('*.counts.tsv')
    )

    script:
    """
    printf 'gene\\tcount\\nENSG0001\\t42\\n' > ${id}.counts.tsv
    """
}
