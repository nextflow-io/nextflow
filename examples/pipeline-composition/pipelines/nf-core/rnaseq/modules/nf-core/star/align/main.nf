nextflow.enable.types = true

process STAR_ALIGN {
    tag "${id}"

    input:
    record(
        id: String,
        fastq_1: Path,
        fastq_2: Path,
        strandedness: String
    )
    fasta: Path
    aligner: String

    output:
    record(
        id: id,
        bam: file('*.bam')
    )

    script:
    """
    echo "aligned ${id} against ${fasta.name} with ${aligner} (${strandedness})" > ${id}.bam
    """
}
