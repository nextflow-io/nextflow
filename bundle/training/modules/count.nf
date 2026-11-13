nextflow.enable.types = true

process COUNT {
    tag id
    label 'small'

    input:
    record(
        id: String,
        group: String,
        part: Path
    )

    output:
    record(
        id: id,
        group: group,
        count: stdout()
    )

    script:
    """
    wc -l < ${part} | tr -d ' \\n'
    """
}
