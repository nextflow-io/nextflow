nextflow.enable.types = true

process SRATOOLS_FASTERQDUMP {
    tag "${id}"

    input:
    id: String

    output:
    record(
        id: id,
        fastq_1: file('*_1.fastq'),
        fastq_2: file('*_2.fastq')
    )

    script:
    """
    echo "@${id}/1" > ${id}_1.fastq
    echo "@${id}/2" > ${id}_2.fastq
    """
}
