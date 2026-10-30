process COUNT {
    tag "${meta.id}"
    label 'small'

    input:
    tuple val(meta), path(part)

    output:
    tuple val(meta), stdout

    script:
    """
    wc -l < ${part} | tr -d ' \\n'
    """
}
